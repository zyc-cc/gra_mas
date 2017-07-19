#ifndef GraphAligner_H
#define GraphAligner_H

//http://biorxiv.org/content/early/2017/04/06/124941
#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/container/flat_set.hpp>
#include <unordered_set>
#include <queue>
#include "vg.pb.h"
#include "SliceRow.h"
#include "SparseBoolMatrix.h"
#include "SparseMatrix.h"
#include "ThreadReadAssertion.h"

using namespace boost;

void printtime(const char* msg)
{
	static auto time = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	auto newtime = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	std::cout << msg << " " << newtime << " (" << (newtime - time) << ")" << std::endl;
	time = newtime;
}

template <typename Word>
class WordConfiguration
{
};

template <>
class WordConfiguration<uint64_t>
{
public:
	static constexpr uint64_t AllZeros = 0x0000000000000000;
	static constexpr uint64_t AllOnes = 0xFFFFFFFFFFFFFFFF;
	static constexpr int WordSize = 64;
	static int popcount(uint64_t x)
	{
		//https://en.wikipedia.org/wiki/Hamming_weight
	    x = (x & 0x5555555555555555) + ((x >>  1) & 0x5555555555555555); //put count of each  2 bits into those  2 bits 
	    x = (x & 0x3333333333333333) + ((x >>  2) & 0x3333333333333333); //put count of each  4 bits into those  4 bits 
	    x = (x & 0x0F0F0F0F0F0F0F0F) + ((x >>  4) & 0x0F0F0F0F0F0F0F0F); //put count of each  8 bits into those  8 bits 
	    x = (x & 0x00FF00FF00FF00FF) + ((x >>  8) & 0x00FF00FF00FF00FF); //put count of each 16 bits into those 16 bits 
	    x = (x & 0x0000FFFF0000FFFF) + ((x >> 16) & 0x0000FFFF0000FFFF); //put count of each 32 bits into those 32 bits 
	    x = (x & 0x00000000FFFFFFFF) + ((x >> 32) & 0x00000000FFFFFFFF); //put count of each 64 bits into those 64 bits 
	    return x;
	}
};

template <typename LengthType, typename ScoreType, typename Word>
class GraphAligner
{
private:
	typedef std::pair<LengthType, LengthType> MatrixPosition;
	class MatrixSlice
	{
	public:
		std::vector<ScoreType> minScorePerWordSlice;
		LengthType finalMinScoreColumn;
		ScoreType finalMinScore;
		size_t cellsProcessed;
	};
	class WordSlice
	{
	public:
		WordSlice() :
		VP(WordConfiguration<Word>::AllZeros),
		VN(WordConfiguration<Word>::AllZeros),
		scoreStart(0),
		scoreEnd(0),
		scoreBeforeStart(0),
		hout(0)
		{}
		WordSlice(Word VP, Word VN, ScoreType scoreStart, ScoreType scoreEnd, ScoreType scoreBeforeStart) :
		VP(VP),
		VN(VN),
		scoreStart(scoreStart),
		scoreEnd(scoreEnd),
		scoreBeforeStart(scoreBeforeStart),
		hout(0)
		{}
		Word VP;
		Word VN;
		ScoreType scoreStart;
		ScoreType scoreEnd;
		ScoreType scoreBeforeStart;
		int hout;
	};
public:
	class SeedHit
	{
	public:
		SeedHit(size_t seqPos, int nodeId, size_t nodePos) : sequencePosition(seqPos), nodeId(nodeId), nodePos(nodePos) {};
		size_t sequencePosition;
		int nodeId;
		size_t nodePos;
	};
	class AlignmentResult
	{
	public:
		AlignmentResult()
		{
		}
		AlignmentResult(vg::Alignment alignment, int maxDistanceFromBand, bool alignmentFailed, size_t cellsProcessed) :
		alignment(alignment),
		maxDistanceFromBand(maxDistanceFromBand),
		alignmentFailed(alignmentFailed),
		cellsProcessed(cellsProcessed)
		{
		}
		vg::Alignment alignment;
		int maxDistanceFromBand;
		bool alignmentFailed;
		size_t cellsProcessed;
	};

	GraphAligner() :
	nodeStart(),
	indexToNode(),
	nodeLookup(),
	nodeIDs(),
	inNeighbors(),
	nodeSequences(),
	gapStartPenalty(1),
	gapContinuePenalty(1),
	finalized(false),
	firstInOrder(0)
	{
		//add the start dummy node as the first node
		dummyNodeStart = nodeSequences.size();
		nodeIDs.push_back(0);
		nodeStart.push_back(nodeSequences.size());
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		reverse.push_back(false);
		nodeSequences.push_back('-');
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
	}
	
	void AddNode(int nodeId, std::string sequence, bool reverseNode)
	{
		assert(!finalized);
		//subgraph extraction might produce different subgraphs with common nodes
		//don't add duplicate nodes
		if (nodeLookup.count(nodeId) != 0) return;

		assert(std::numeric_limits<LengthType>::max() - sequence.size() > nodeSequences.size());
		nodeLookup[nodeId] = nodeStart.size();
		nodeIDs.push_back(nodeId);
		nodeStart.push_back(nodeSequences.size());
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		reverse.push_back(reverseNode);
		nodeSequences.insert(nodeSequences.end(), sequence.begin(), sequence.end());
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		assert(nodeIDs.size() == nodeStart.size());
		assert(nodeStart.size() == inNeighbors.size());
		assert(inNeighbors.size() == nodeEnd.size());
		assert(nodeEnd.size() == notInOrder.size());
		assert(nodeSequences.size() == indexToNode.size());
		assert(inNeighbors.size() == outNeighbors.size());
	}
	
	void AddEdgeNodeId(int node_id_from, int node_id_to)
	{
		assert(!finalized);
		assert(nodeLookup.count(node_id_from) > 0);
		assert(nodeLookup.count(node_id_to) > 0);
		auto from = nodeLookup[node_id_from];
		auto to = nodeLookup[node_id_to];
		assert(to >= 0);
		assert(from >= 0);
		assert(to < inNeighbors.size());
		assert(from < nodeStart.size());

		//subgraph extraction might produce different subgraphs with common edges
		//don't add duplicate edges
		if (std::find(inNeighbors[to].begin(), inNeighbors[to].end(), from) != inNeighbors[to].end()) return;

		inNeighbors[to].push_back(from);
		outNeighbors[from].push_back(to);
		if (from >= to)
		{
			notInOrder[to] = true;
		}
	}

	void Finalize()
	{
		//add the end dummy node as the last node
		dummyNodeEnd = nodeSequences.size();
		nodeIDs.push_back(0);
		nodeStart.push_back(nodeSequences.size());
		reverse.push_back(false);
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		nodeSequences.push_back('-');
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		assert(nodeSequences.size() >= nodeStart.size());
		assert(nodeEnd.size() == nodeStart.size());
		assert(notInOrder.size() == nodeStart.size());
		assert(inNeighbors.size() == nodeStart.size());
		assert(outNeighbors.size() == nodeStart.size());
		assert(notInOrder.size() == nodeStart.size());
		assert(reverse.size() == nodeStart.size());
		assert(nodeIDs.size() == nodeStart.size());
		assert(indexToNode.size() == nodeSequences.size());
		finalized = true;
		int specialNodes = 0;
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			if (inNeighbors[i].size() >= 2) specialNodes++;
		}
		std::cerr << specialNodes << " nodes with in-degree >= 2" << std::endl;
		firstInOrder = 0;
		for (size_t i = 1; i < notInOrder.size(); i++)
		{
			if (notInOrder[i]) firstInOrder = i+1;
			//all not-in-order nodes have to be at the start
			assert(i == 1 || !notInOrder[i] || notInOrder[i-1]);
		}
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart) const
	{
		assert(finalized);
		auto band = getFullBand(sequence.size(), dynamicRowStart);
		auto trace = getBacktrace(sequence, dynamicWidth, dynamicRowStart, band);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace), std::get<3>(trace));
		return result;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart, const std::vector<SeedHit>& seedHits, int startBandwidth) const
	{
		assert(finalized);
		auto band = getSeededStartBand(seedHits, dynamicRowStart, startBandwidth, sequence);
		auto trace = getBacktrace(sequence, dynamicWidth, dynamicRowStart, band);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace), std::get<3>(trace));
		return result;
	}

	size_t SizeInBp()
	{
		return nodeSequences.size();
	}

private:

	class ExpandoCell
	{
	public:
		ExpandoCell(LengthType w, LengthType j, size_t bt) :
		position(w, j),
		backtraceIndex(bt)
		{}
		MatrixPosition position;
		size_t backtraceIndex;
	};

	std::vector<MatrixPosition> backtrace(MatrixPosition endPosition, const std::string& sequence, const std::vector<ScoreType>& minScorePerWordSlice) const
	{
		//first try the faster way
		auto result = backtraceInner(endPosition, sequence, minScorePerWordSlice, false);
		//if that breaks then do the full backtrace
		if (result.size() == 0) result = backtraceInner(endPosition, sequence, minScorePerWordSlice, true);
		assert(result.size() != 0);
		return result;
	}

	std::vector<MatrixPosition> backtraceInner(MatrixPosition endPosition, const std::string& sequence, const std::vector<ScoreType>& minScorePerWordSlice, bool fullBacktrace) const
	{
		ScoreType scoreAtEnd = minScorePerWordSlice.back()+1;
		ScoreType currentDistance = 0;
		std::vector<ExpandoCell> visitedExpandos;
		std::vector<ExpandoCell> currentDistanceQueue;
		std::vector<ExpandoCell> currentDistancePlusOneQueue;
		currentDistanceQueue.emplace_back(endPosition.first, endPosition.second, 0);
		SparseBoolMatrix<SliceRow<LengthType>> visitedCells {nodeSequences.size(), sequence.size()+1};
		while (true)
		{
			if (currentDistanceQueue.size() == 0)
			{
				if (currentDistancePlusOneQueue.size() == 0)
				{
					return std::vector<MatrixPosition> {};
				}
				assert(currentDistancePlusOneQueue.size() > 0);
				std::swap(currentDistanceQueue, currentDistancePlusOneQueue);
				currentDistance++;
				if (currentDistance > scoreAtEnd)
				{
					return std::vector<MatrixPosition> {};
				}
			}
			auto current = currentDistanceQueue.back();
			currentDistanceQueue.pop_back();
			auto w = current.position.first;
			auto j = current.position.second;
			if (j == 0)
			{
				visitedExpandos.push_back(current);
				break;
			}
			auto sliceIndex = (j-1) / WordConfiguration<Word>::WordSize;
			assert(sliceIndex < minScorePerWordSlice.size());
			ScoreType maxDistanceHere;
			if (fullBacktrace)
			{
				maxDistanceHere = scoreAtEnd;
			}
			else
			{
				maxDistanceHere = scoreAtEnd - minScorePerWordSlice[sliceIndex];
			}
			if (currentDistance > maxDistanceHere) continue;
			if (visitedCells.get(w, j)) continue;
			visitedCells.set(w, j);
			visitedExpandos.push_back(current);
			auto nodeIndex = indexToNode[w];
			auto backtraceIndexToCurrent = visitedExpandos.size()-1;
			currentDistancePlusOneQueue.emplace_back(w, j-1, backtraceIndexToCurrent);
			if (w == nodeStart[nodeIndex])
			{
				for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
				{
					auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
					currentDistancePlusOneQueue.emplace_back(u, j, backtraceIndexToCurrent);
					if (sequence[j-1] == 'N' || nodeSequences[w] == sequence[j-1])
					{
						currentDistanceQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
					}
					else
					{
						currentDistancePlusOneQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
					}
				}
			}
			else
			{
				auto u = w-1;
				currentDistancePlusOneQueue.emplace_back(u, j, backtraceIndexToCurrent);
				if (sequence[j-1] == 'N' || nodeSequences[w] == sequence[j-1])
				{
					currentDistanceQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
				}
				else
				{
					currentDistancePlusOneQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
				}
			}
		}
		std::cerr << "backtrace visited " << visitedCells.totalOnes() << " cells" << std::endl;
		assert(currentDistance <= scoreAtEnd);
		auto index = visitedExpandos.size()-1;
		std::vector<MatrixPosition> result;
		while (index > 0)
		{
			result.push_back(visitedExpandos[index].position);
			assert(visitedExpandos[index].backtraceIndex < index);
			index = visitedExpandos[index].backtraceIndex;
		}
		return result;
	}

	std::vector<std::vector<bool>> getFullBand(size_t sequenceSize, LengthType dynamicRowStart) const
	{
		std::vector<std::vector<bool>> result;
		result.resize(dynamicRowStart/WordConfiguration<Word>::WordSize);
		for (size_t i = 0; i < dynamicRowStart/WordConfiguration<Word>::WordSize; i++)
		{
			result[i].resize(nodeStart.size(), true);
		}
		return result;
	}

	std::vector<std::vector<bool>> getSeededStartBand(const std::vector<SeedHit>& originalSeedHits, int dynamicRowStart, int startBandwidth, const std::string& sequence) const
	{
		auto seedHitsInMatrix = getSeedHitPositionsInMatrix(sequence, originalSeedHits);
		auto bandLocations = getBandLocations(sequence.size(), seedHitsInMatrix, dynamicRowStart);
		std::vector<std::vector<bool>> result;
		result.resize(dynamicRowStart/WordConfiguration<Word>::WordSize);
		for (LengthType j = 0; j < dynamicRowStart && j < sequence.size()+1; j += WordConfiguration<Word>::WordSize)
		{
			auto index = j/WordConfiguration<Word>::WordSize;
			result[index].resize(nodeStart.size(), false);
			for (size_t i = 0; i < bandLocations[j].size(); i++)
			{
				expandBandDynamically(result[index], bandLocations[j][i], startBandwidth);
			}
		}
		return result;
	}

	std::vector<std::vector<LengthType>> getBandLocations(int sequenceLength, const std::vector<MatrixPosition>& seedHits, LengthType maxRow) const
	{
		if (maxRow > sequenceLength+1) maxRow = sequenceLength+1;
		std::vector<std::vector<LengthType>> forwardResult;
		std::vector<std::vector<LengthType>> backwardResult;
		backwardResult.resize(sequenceLength+1);
		forwardResult.resize(sequenceLength+1);
		backwardResult[0].emplace_back(0);
		forwardResult[0].emplace_back(0);
		for (auto hit : seedHits)
		{
			if (hit.second < maxRow) expandBandForwards(forwardResult, hit.first, hit.second, sequenceLength, maxRow);
			expandBandBackwards(backwardResult, hit.first, hit.second, sequenceLength);
		}
		std::vector<std::vector<LengthType>> result;
		result.resize(maxRow+1);
		for (size_t j = 0; j < maxRow; j++)
		{
			std::set<LengthType> rowResult;
			rowResult.insert(forwardResult[j].begin(), forwardResult[j].end());
			rowResult.insert(backwardResult[j].begin(), backwardResult[j].end());
			result[j].insert(result[j].end(), rowResult.begin(), rowResult.end());
		}
		return result;
	}

	std::vector<MatrixPosition> getSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const
	{
		std::vector<MatrixPosition> result;
		for (size_t i = 0; i < seedHits.size(); i++)
		{
			assert(nodeLookup.count(seedHits[i].nodeId) > 0);
			result.emplace_back(nodeStart[nodeLookup.at(seedHits[i].nodeId)] + seedHits[i].nodePos, seedHits[i].sequencePosition);
		}
		return result;
	}

	AlignmentResult emptyAlignment() const
	{
		vg::Alignment result;
		result.set_score(std::numeric_limits<decltype(result.score())>::min());
		return AlignmentResult { result, 0, true, 0 };
	}

	AlignmentResult traceToAlignment(const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, int maxDistanceFromBand, size_t cellsProcessed) const
	{
		vg::Alignment result;
		result.set_name(seq_id);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		size_t pos = 0;
		size_t oldNode = indexToNode[trace[0].first];
		while (oldNode == dummyNodeStart)
		{
			pos++;
			if (pos == trace.size()) return emptyAlignment();
			assert(pos < trace.size());
			oldNode = indexToNode[trace[pos].first];
			assert(oldNode < nodeIDs.size());
		}
		if (oldNode == dummyNodeEnd) return emptyAlignment();
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		position->set_node_id(nodeIDs[oldNode]);
		position->set_is_reverse(reverse[oldNode]);
		for (; pos < trace.size(); pos++)
		{
			if (indexToNode[trace[pos].first] == dummyNodeEnd) break;
			if (indexToNode[trace[pos].first] == oldNode) continue;
			oldNode = indexToNode[trace[pos].first];
			rank++;
			vgmapping = path->add_mapping();
			position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(rank);
			position->set_node_id(nodeIDs[oldNode]);
			position->set_is_reverse(reverse[oldNode]);
		}
		result.set_score(score);
		result.set_sequence(sequence);
		return AlignmentResult { result, maxDistanceFromBand, false, cellsProcessed };
	}

	void expandBandForwards(std::vector<std::vector<LengthType>>& result, LengthType w, LengthType j, size_t sequenceLength, LengthType maxRow) const
	{
		if (std::find(result[j].begin(), result[j].end(), w) != result[j].end()) return;
		auto nodeIndex = indexToNode[w];
		auto end = nodeEnd[nodeIndex];
		while (w != end && j < sequenceLength+1 && j < maxRow)
		{
			result[j].emplace_back(w);
			w++;
			j++;
		}
		if (w == end && j < sequenceLength+1 && j < maxRow)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandBandForwards(result, nodeStart[outNeighbors[nodeIndex][i]], j, sequenceLength, maxRow);
			}
		}
	}

	void expandBandBackwards(std::vector<std::vector<LengthType>>& result, LengthType w, LengthType j, size_t sequenceLength) const
	{
		if (std::find(result[j].begin(), result[j].end(), w) != result[j].end()) return;
		auto nodeIndex = indexToNode[w];
		auto start = nodeStart[nodeIndex];
		while (w != start && j > 0)
		{
			result[j].emplace_back(w);
			w--;
			j--;
		}
		result[j].emplace_back(w);
		if (w == start && j > 0)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandBandBackwards(result, nodeEnd[inNeighbors[nodeIndex][i]] - 1, j-1, sequenceLength);
			}
		}
	}

	class IndexWithScore
	{
	public:
		IndexWithScore(LengthType index, ScoreType score) :
		index(index),
		score(score)
		{}
		LengthType index;
		ScoreType score;
		bool operator>(const IndexWithScore& other) const
		{
			return score > other.score;
		}
	};

	void expandBandFromPrevious(std::vector<bool>& currentBand, const std::vector<bool>& previousBand, const std::vector<WordSlice>& previousSlice, LengthType dynamicWidth, const std::vector<ScoreType>& nodeMinScores) const
	{
		assert(dynamicWidth < nodeSequences.size());
		assert(currentBand.size() == previousBand.size());
		assert(currentBand.size() == nodeStart.size());
		std::priority_queue<IndexWithScore, std::vector<IndexWithScore>, std::greater<IndexWithScore>> nodeQueue;
		std::priority_queue<IndexWithScore, std::vector<IndexWithScore>, std::greater<IndexWithScore>> endQueue;
		for (size_t i = 0; i < nodeStart.size(); i++)
		{
			if (!previousBand[i]) continue;
			nodeQueue.emplace(i, nodeMinScores[i]);
			endQueue.emplace(i, previousSlice[nodeEnd[i]-1].scoreEnd);
		}
		LengthType currentWidth = 0;
		while (currentWidth < dynamicWidth)
		{
			assert(nodeQueue.size() > 0 || endQueue.size() > 0);
			if (nodeQueue.size() == 0)
			{
				auto nextNode = endQueue.top();
				endQueue.pop();
				for (size_t i = 0; i < outNeighbors[nextNode.index].size(); i++)
				{
					auto neighbor = outNeighbors[nextNode.index][i];
					if (!currentBand[neighbor])
					{
						currentBand[neighbor] = true;
						currentWidth += nodeEnd[neighbor] - nodeStart[neighbor];
						endQueue.emplace(neighbor, nextNode.score + nodeEnd[neighbor] - nodeStart[neighbor]);
					}
				}
				continue;
			}
			if (endQueue.size() == 0)
			{
				auto nextNode = nodeQueue.top();
				nodeQueue.pop();
				if (!currentBand[nextNode.index])
				{
					currentBand[nextNode.index] = true;
					currentWidth += nodeEnd[nextNode.index] - nodeStart[nextNode.index];
				}
				continue;
			}
			auto nodeBest = nodeQueue.top();
			auto endBest = endQueue.top();
			if (nodeBest.score <= endBest.score)
			{
				nodeQueue.pop();
				assert(previousBand[nodeBest.index]);
				if (!currentBand[nodeBest.index])
				{
					currentBand[nodeBest.index] = true;
					currentWidth += nodeEnd[nodeBest.index] - nodeStart[nodeBest.index];
				}
			}
			else
			{
				endQueue.pop();
				assert(currentBand[endBest.index]);
				for (size_t i = 0; i < outNeighbors[endBest.index].size(); i++)
				{
					auto neighbor = outNeighbors[endBest.index][i];
					if (!currentBand[neighbor])
					{
						currentBand[neighbor] = true;
						currentWidth += nodeEnd[neighbor] - nodeStart[neighbor];
						endQueue.emplace(neighbor, endBest.score + nodeEnd[neighbor] - nodeStart[neighbor]);
					}
				}
			}
		}
	}

	void expandBandDynamically(std::vector<bool>& band, LengthType previousMinimumIndex, LengthType dynamicWidth) const
	{
		assert(previousMinimumIndex < nodeSequences.size());
		auto nodeIndex = indexToNode[previousMinimumIndex];
		band[nodeIndex] = true;
		LengthType start = nodeStart[nodeIndex];
		LengthType end = nodeEnd[nodeIndex];
		if (dynamicWidth > previousMinimumIndex - start)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], dynamicWidth - (previousMinimumIndex - start));
			}
		}
		if (dynamicWidth > end - previousMinimumIndex)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandForward(band, outNeighbors[nodeIndex][i], dynamicWidth - (end - previousMinimumIndex));
			}
		}
	}

	void expandDynamicBandBackward(std::vector<bool>& band, LengthType nodeIndex, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < nodeStart.size());
		//todo fix: currently only the first path that reaches the node is considered
		//this means that it might not reach some nodes with distance < dynamicwidth
		if (band[nodeIndex]) return;
		band[nodeIndex] = true;
		for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
		{
			expandDynamicBandForward(band, outNeighbors[nodeIndex][i], dynamicWidth - 1);
		}
		auto nodeSize = nodeEnd[nodeIndex] - nodeStart[nodeIndex];
		if (dynamicWidth > nodeSize)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], dynamicWidth - nodeSize);
			}
		}
	}

	template <typename MatrixType>
	void expandDynamicBandForward(MatrixType& band, LengthType nodeIndex, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < nodeStart.size());
		//todo fix: currently only the first path that reaches the node is considered
		//this means that it might not reach some nodes with distance < dynamicwidth if there's multiple paths and it arbitrarily picks the longest one first
		if (band[nodeIndex]) return;
		band[nodeIndex] = true;
		for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
		{
			expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], dynamicWidth - 1);
		}
		auto nodeSize = nodeEnd[nodeIndex] - nodeStart[nodeIndex];
		if (dynamicWidth > nodeSize)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandForward(band, outNeighbors[nodeIndex][i], dynamicWidth - nodeSize);
			}
		}
	}

	uint64_t bytePopcounts(uint64_t value) const
	{
		uint64_t x = value;
		x -= (x >> 1) & 0x5555555555555555;
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
		return x;
	}

	uint64_t bytePrefixSums(uint64_t value, int addition) const
	{
		value <<= 8;
		assert(addition >= 0);
		value += addition;
		return value * 0x0101010101010101;
	}

	uint64_t byteVPVNSum(uint64_t prefixSumVP, uint64_t prefixSumVN) const
	{
		uint64_t result = 0x8080808080808080;
		assert((prefixSumVP & result) == 0);
		assert((prefixSumVN & result) == 0);
		result += prefixSumVP;
		result -= prefixSumVN;
		return result;
	}

	std::string wordToStr(uint64_t value) const
	{
		std::string result;
		for (int i = 63; i >= 0; i--)
		{
			if (value & (((uint64_t)1) << i))
			{
				result += "1";
			}
			else
			{
				result += "0";
			}
		}
		return result;
	}

	std::pair<uint64_t, uint64_t> differenceMasks(uint64_t leftVP, uint64_t leftVN, uint64_t rightVP, uint64_t rightVN, int scoreDifference) const
	{
		//do the addition in chunks to prevent calculations from interfering with the other bytes
		const uint64_t chunkmask1 = 0xFF00FF00FF00FF00;
		const uint64_t chunkmask2 = 0x00FF00FF00FF00FF;
		//the fences make sure that when moving positive<->negative the other bytes are not effected
		const uint64_t fencemask1 = 0x0055005500550055;
		const uint64_t fencemask2 = 0x5500550055005500;
		const uint64_t signmask = 0x8080808080808080;
		assert(scoreDifference >= 0);
		uint64_t byteVPVNSumLeft = byteVPVNSum(bytePrefixSums(bytePopcounts(leftVP), 0), bytePrefixSums(bytePopcounts(leftVN), 0));
		uint64_t byteVPVNSumRight = byteVPVNSum(bytePrefixSums(bytePopcounts(rightVP), scoreDifference), bytePrefixSums(bytePopcounts(rightVN), 0));
		uint64_t left = (byteVPVNSumLeft & chunkmask1) | fencemask1;
		uint64_t right = (byteVPVNSumRight & chunkmask1);
		uint64_t difference = (left - right) & chunkmask1;
		left = (byteVPVNSumLeft & chunkmask2) | fencemask2;
		right = (byteVPVNSumRight & chunkmask2);
		difference |= (left - right) & chunkmask2;
		//difference now contains the prefix sum difference (left-right) at each byte
		uint64_t resultLeftSmallerThanRight = 0;
		uint64_t resultRightSmallerThanLeft = 0;
		for (int bit = 0; bit < 8; bit++)
		{
			uint64_t bitmask = 0x0101010101010101 << bit;
			uint64_t difference1 = (difference & chunkmask1) | fencemask1;
			uint64_t difference2 = (difference & chunkmask2) | fencemask2;
			difference1 += (leftVP & bitmask & chunkmask1) >> bit;
			difference1 -= (leftVN & bitmask & chunkmask1) >> bit;
			difference1 -= (rightVP & bitmask & chunkmask1) >> bit;
			difference1 += (rightVN & bitmask & chunkmask1) >> bit;
			difference2 += (leftVP & bitmask & chunkmask2) >> bit;
			difference2 -= (leftVN & bitmask & chunkmask2) >> bit;
			difference2 -= (rightVP & bitmask & chunkmask2) >> bit;
			difference2 += (rightVN & bitmask & chunkmask2) >> bit;
			difference = (difference1 & chunkmask1) | (difference2 & chunkmask2);
			//difference now contains the prefix sums difference (left-right) at each byte at (bit)'th bit
			//left < right when the prefix sum difference is negative (sign bit is set)
			uint64_t negative = (difference & signmask);
			resultLeftSmallerThanRight |= negative >> (7 - bit);
			//Test equality to zero. If it's zero, substracting one will make the sign bit 0, otherwise 1
			uint64_t notEqualToZero = ((difference | signmask) - 0x0101010101010101) & signmask;
			//right > left when the prefix sum difference is positive (not zero and not negative)
			resultRightSmallerThanLeft |= (notEqualToZero & ~negative) >> (7 - bit);
		}
		return std::make_pair(resultLeftSmallerThanRight, resultRightSmallerThanLeft);
	}

	WordSlice mergeTwoSlices(WordSlice left, WordSlice right) const
	{
		//O(log w), because prefix sums need log w chunks of log w bits
		static_assert(std::is_same<Word, uint64_t>::value);
		if (left.scoreBeforeStart > right.scoreBeforeStart) std::swap(left, right);
		WordSlice result;
		assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
		assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
		// std::cerr << "left.scoreBeforeStart " << left.scoreBeforeStart << std::endl;
		// std::cerr << "right.scoreBeforeStart " << right.scoreBeforeStart << std::endl;
		auto masks = differenceMasks(left.VP, left.VN, right.VP, right.VN, right.scoreBeforeStart - left.scoreBeforeStart);
		auto leftSmaller = masks.first;
		auto rightSmaller = masks.second;
		auto mask = (rightSmaller | ((leftSmaller | rightSmaller) - (rightSmaller << 1))) & ~leftSmaller;
		// std::cerr << "left.VP          " << wordToStr(left.VP) << std::endl;
		// std::cerr << "left.VN          " << wordToStr(left.VN) << std::endl;
		// std::cerr << "right.VP         " << wordToStr(right.VP) << std::endl;
		// std::cerr << "right.VN         " << wordToStr(right.VN) << std::endl;
		// std::cerr << "leftSmaller      " << wordToStr(leftSmaller) << std::endl;
		// std::cerr << "rightSmaller     " << wordToStr(rightSmaller) << std::endl;
		// std::cerr << "mask             " << wordToStr(mask) << std::endl;
		uint64_t leftLow = ~left.VN;
		uint64_t leftHigh = left.VP;
		uint64_t rightLow = ~right.VN;
		uint64_t rightHigh = right.VP;
		// std::cerr << "leftHigh(pre)    " << wordToStr(leftHigh) << std::endl;
		// std::cerr << "leftLow(pre)     " << wordToStr(leftLow) << std::endl;
		// std::cerr << "rightHigh(pre)   " << wordToStr(rightHigh) << std::endl;
		// std::cerr << "rightLow(pre)    " << wordToStr(rightLow) << std::endl;
		uint64_t leftReduction = leftSmaller & (rightSmaller << 1);
		uint64_t rightReduction = rightSmaller & (leftSmaller << 1);
		if ((rightSmaller & 1) && left.scoreBeforeStart < right.scoreBeforeStart)
		{
			rightReduction |= 1;
		}
		// std::cerr << "leftReduction    " << wordToStr(leftReduction) << std::endl;
		// std::cerr << "rightReduction   " << wordToStr(rightReduction) << std::endl;
		assert((leftReduction & rightHigh) == leftReduction);
		assert((rightReduction & leftHigh) == rightReduction);
		assert((leftReduction & leftLow) == 0);
		assert((rightReduction & rightLow) == 0);
		leftLow |= leftReduction;
		rightLow |= rightReduction;
		// std::cerr << "leftHigh(post)   " << wordToStr(leftHigh) << std::endl;
		// std::cerr << "leftLow(post)    " << wordToStr(leftLow) << std::endl;
		// std::cerr << "rightHigh(post)  " << wordToStr(rightHigh) << std::endl;
		// std::cerr << "rightLow(post)   " << wordToStr(rightLow) << std::endl;
		uint64_t resultLow = (leftLow & ~mask) | (rightLow & mask);
		uint64_t resultHigh = (leftHigh & ~mask) | (rightHigh & mask);
		// std::cerr << "resultHigh       " << wordToStr(resultHigh) << std::endl;
		// std::cerr << "resultLow        " << wordToStr(resultLow) << std::endl;
		assert((resultHigh & ~resultLow) == 0);
		result.VP = resultHigh;
		result.VN = ~resultLow;
		// std::cerr << "resultVP         " << wordToStr(result.VP) << std::endl;
		// std::cerr << "resultVN         " << wordToStr(result.VN) << std::endl;
		result.scoreBeforeStart = std::min(left.scoreBeforeStart, right.scoreBeforeStart);
		result.scoreStart = std::min(left.scoreStart, right.scoreStart);
		result.scoreEnd = std::min(left.scoreEnd, right.scoreEnd);
#ifndef NDEBUG
		auto correctValue = mergeTwoSlicesCellByCell(left, right);
		// std::cerr << "correctVP        " << wordToStr(correctValue.VP) << std::endl;
		// std::cerr << "correctVN        " << wordToStr(correctValue.VN) << std::endl;
		assert(result.VP == correctValue.VP);
		assert(result.VN == correctValue.VN);
		assert(result.scoreBeforeStart == correctValue.scoreBeforeStart);
		assert(result.scoreStart == correctValue.scoreStart);
		assert(result.scoreEnd == correctValue.scoreEnd);
#endif
		return result;
	}

	WordSlice mergeTwoSlicesCellByCell(WordSlice left, WordSlice right) const
	{
		//currently O(w), O(log^2 w) is possible
		//todo figure out the details and implement
		//todo is O(log w) possible?
		assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
		assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
		ScoreType leftScore = left.scoreStart;
		WordSlice merged;
		merged.scoreBeforeStart = std::min(left.scoreBeforeStart, right.scoreBeforeStart);
		merged.VP = WordConfiguration<Word>::AllZeros;
		merged.VN = WordConfiguration<Word>::AllZeros;
		ScoreType rightScore = right.scoreStart;
		merged.scoreStart = std::min(rightScore, leftScore);
		assert(merged.scoreStart >= merged.scoreBeforeStart - 1);
		assert(merged.scoreStart <= merged.scoreBeforeStart + 1);
		if (merged.scoreStart == merged.scoreBeforeStart - 1)
		{
			merged.VN |= 1;
		}
		else if (merged.scoreStart == merged.scoreBeforeStart + 1)
		{
			merged.VP |= 1;
		}
		ScoreType previousScore = merged.scoreStart;
		for (size_t j = 1; j < WordConfiguration<Word>::WordSize; j++)
		{
			Word mask = ((Word)1) << j;
			if (left.VP & mask) leftScore++;
			else if (left.VN & mask) leftScore--;
			if (right.VN & mask) rightScore--;
			else if (right.VP & mask) rightScore++;
			ScoreType betterScore = std::min(leftScore, rightScore);
			if (betterScore == previousScore+1) merged.VP |= mask;
			else if (betterScore == previousScore-1) merged.VN |= mask;
			assert((merged.VP & merged.VN) == WordConfiguration<Word>::AllZeros);
			assert(betterScore >= previousScore-1);
			assert(betterScore <= previousScore+1);
			previousScore = betterScore;
		}
		merged.scoreEnd = previousScore;
		assert((merged.VP & merged.VN) == WordConfiguration<Word>::AllZeros);
		assert(merged.scoreEnd <= left.scoreEnd);
		assert(merged.scoreEnd <= right.scoreEnd);
		assert(merged.scoreStart <= left.scoreStart);
		assert(merged.scoreStart <= right.scoreStart);
		return merged;
	}

	WordSlice getNodeStartSlice(Word Eq, size_t nodeIndex, const std::vector<WordSlice>& previousSlice, const std::vector<WordSlice>& currentSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		WordSlice previous;
		bool foundOne = false;
		bool forceFirstHorizontalPositive = false;
		int hin = 0;
		if (previousBand[nodeIndex]) hin = previousSlice[nodeStart[nodeIndex]].hout;
		for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
		{
			auto neighbor = inNeighbors[nodeIndex][i];
			if (!currentBand[neighbor]) continue;
			if (!foundOne)
			{
				previous = currentSlice[nodeEnd[neighbor]-1];
				forceFirstHorizontalPositive = firstZeroForced(previousBand, neighbor, nodeEnd[neighbor]-1, Eq, currentSlice, previousSlice);
				foundOne = true;
			}
			else
			{
				auto competitor = currentSlice[nodeEnd[neighbor]-1];
				if (competitor.scoreBeforeStart < previous.scoreBeforeStart)
				{
					forceFirstHorizontalPositive = firstZeroForced(previousBand, neighbor, nodeEnd[neighbor]-1, Eq, currentSlice, previousSlice);
				}
				else if (competitor.scoreBeforeStart == previous.scoreBeforeStart)
				{
					forceFirstHorizontalPositive = forceFirstHorizontalPositive & firstZeroForced(previousBand, neighbor, nodeEnd[neighbor]-1, Eq, currentSlice, previousSlice);
				}
				previous = mergeTwoSlices(previous, competitor);
			}
		}
		assert(foundOne);
		auto result = getNextSlice(Eq, previous, hin, forceFirstHorizontalPositive);
		return result;
	}

	WordSlice getSourceSliceWithoutBefore(size_t row) const
	{
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, row+1, row+WordConfiguration<Word>::WordSize, row };
	}

	WordSlice getSourceSliceFromColumn(size_t column, const std::vector<WordSlice>& previousSlice) const
	{
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousSlice[column].scoreEnd+1, previousSlice[column].scoreEnd+WordConfiguration<Word>::WordSize, previousSlice[column].scoreEnd };
	}

	WordSlice getSourceSlice(size_t nodeIndex, const std::vector<WordSlice>& previousSlice) const
	{
		auto start = nodeStart[nodeIndex];
		return getSourceSliceFromColumn(start, previousSlice);
	}

	bool isSource(size_t nodeIndex, const std::vector<bool>& band) const
	{
		for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
		{
			if (band[inNeighbors[nodeIndex][i]]) return false;
		}
		return true;
	}

	Word getEq(Word BA, Word BT, Word BC, Word BG, LengthType w) const
	{
		switch(nodeSequences[w])
		{
			case 'A':
				return BA;
				break;
			case 'T':
				return BT;
				break;
			case 'C':
				return BC;
				break;
			case 'G':
				return BG;
				break;
			case '-':
				assert(false);
				break;
			default:
				assert(false);
				break;
		}
		assert(false);
		return 0;
	}

	WordSlice getNextSlice(Word Eq, WordSlice slice, int hin, bool forceFirstHorizontalPositive) const
	{
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408
		Word Xv = Eq | slice.VN;
		//between 7 and 8
		if (hin < 0) Eq |= 1;
		if (forceFirstHorizontalPositive)
		{
			Eq &= ~((Word)1);
		}
		Word Xh = (((Eq & slice.VP) + slice.VP) ^ slice.VP) | Eq;
		Word Ph = slice.VN | ~(Xh | slice.VP);
		Word Mh = slice.VP & Xh;
		if (forceFirstHorizontalPositive)
		{
			Ph |= 1;
			Mh &= ~((Word)1);
		}
		if (Ph & ((Word)1))
		{
			slice.scoreStart += 1;
		}
		else if (Mh & ((Word)1))
		{
			slice.scoreStart -= 1;
		}
		Word lastBitMask = (((Word)1) << (WordConfiguration<Word>::WordSize - 1));
		if (Ph & lastBitMask)
		{
			slice.scoreEnd += 1;
			slice.hout = 1;
		}
		else if (Mh & lastBitMask)
		{
			slice.scoreEnd -= 1;
			slice.hout = -1;
		}
		else
		{
			slice.hout = 0;
		}
		Ph <<= 1;
		Mh <<= 1;
		//between 16 and 17
		if (hin < 0) Mh |= 1; else if (hin > 0) Ph |= 1;
		slice.VP = Mh | ~(Xv | Ph);
		slice.VN = Ph & Xv;

#ifndef NDEBUG
		auto wcvpExceptFirst = WordConfiguration<Word>::popcount(slice.VP & ~((Word)1));
		auto wcvnExceptFirst = WordConfiguration<Word>::popcount(slice.VN & ~((Word)1));
		assert(slice.scoreEnd == slice.scoreStart + wcvpExceptFirst - wcvnExceptFirst);
#endif
		slice.scoreBeforeStart = slice.scoreStart;
		if (slice.VP & 1) slice.scoreBeforeStart -= 1;
		if (slice.VN & 1) slice.scoreBeforeStart += 1;

		return slice;
	}

	class NodeCalculationResult
	{
	public:
		ScoreType minScore;
		LengthType minScoreIndex;
		size_t cellsProcessed;
	};

	bool firstZeroForced(const std::vector<bool>& previousBand, LengthType neighborNodeIndex, LengthType neighborColumn, Word currentEq, const std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice) const
	{
		if (previousBand[neighborNodeIndex])
		{
			if (currentSlice[neighborColumn].scoreStart < previousSlice[neighborColumn].scoreEnd)
			{
				return true;
			}
			if (currentSlice[neighborColumn].scoreStart == previousSlice[neighborColumn].scoreEnd && !(currentEq & 1))
			{
				return true;
			}
			return false;
		}
		else
		{
			return true;
		}
		assert(false);
	}

	NodeCalculationResult calculateNode(size_t i, size_t j, Word BA, Word BT, Word BC, Word BG, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, bool forceSource) const
	{
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreIndex = 0;
		result.cellsProcessed = 0;

		LengthType start = nodeStart[i];
		if (forceSource || isSource(i, currentBand))
		{
			if (previousBand[i])
			{
				currentSlice[start] = getSourceSlice(i, previousSlice);
			}
			else
			{
				currentSlice[start] = getSourceSliceWithoutBefore(j);
			}
			if (currentSlice[start].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[start].scoreEnd;
				result.minScoreIndex = start;
			}
			start++;
		}
		else
		{
			Word Eq = getEq(BA, BT, BC, BG, nodeStart[i]);
			currentSlice[start] = getNodeStartSlice(Eq, i, previousSlice, currentSlice, currentBand, previousBand);
			if (previousBand[i] && currentSlice[start].scoreBeforeStart > previousSlice[start].scoreEnd)
			{
				currentSlice[start] = mergeTwoSlices(getSourceSliceFromColumn(start, previousSlice), currentSlice[start]);
			}
			if (currentSlice[start].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[start].scoreEnd;
				result.minScoreIndex = start;
			}
			start++;
		}

		assert(start == nodeStart[i]+1);

		for (LengthType w = start; w < nodeEnd[i]; w++)
		{
			Word Eq = getEq(BA, BT, BC, BG, w);
			int hin = 0;
			if (previousBand[i]) hin = previousSlice[w].hout;
			bool forceFirstHorizontalPositive = firstZeroForced(previousBand, i, w-1, Eq, currentSlice, previousSlice);

			currentSlice[w] = getNextSlice(Eq, currentSlice[w-1], hin, forceFirstHorizontalPositive);

			if (previousBand[i] && currentSlice[w].scoreBeforeStart > previousSlice[w].scoreEnd)
			{
				currentSlice[w] = mergeTwoSlices(getSourceSliceFromColumn(w, previousSlice), currentSlice[w]);
			}
			if (currentSlice[w].scoreBeforeStart > j)
			{
				currentSlice[w] = mergeTwoSlices(getSourceSliceWithoutBefore(j), currentSlice[w]);
			}

#ifndef NDEBUG
			auto wcvp = WordConfiguration<Word>::popcount(currentSlice[w].VP);
			auto wcvn = WordConfiguration<Word>::popcount(currentSlice[w].VN);
			assert(currentSlice[w].scoreEnd == currentSlice[w].scoreBeforeStart + wcvp - wcvn);
#endif
			assert(previousBand[i] || currentSlice[w].scoreBeforeStart == j || currentSlice[w].scoreStart == currentSlice[w-1].scoreStart + 1);
			assert(currentSlice[w].scoreStart >= 0);
			assert(currentSlice[w].scoreEnd >= 0);
			assert(currentSlice[w].scoreStart <= currentSlice[w].scoreEnd + WordConfiguration<Word>::WordSize);
			assert(currentSlice[w].scoreEnd <= currentSlice[w].scoreStart + WordConfiguration<Word>::WordSize);
			assert((currentSlice[w].VP & currentSlice[w].VN) == WordConfiguration<Word>::AllZeros);

			assert(!previousBand[i] || currentSlice[w].scoreBeforeStart <= previousSlice[w].scoreEnd);
			assert(currentSlice[w].scoreBeforeStart >= 0);
			if (currentSlice[w].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[w].scoreEnd;
				result.minScoreIndex = w;
			}
		}
		result.cellsProcessed = (nodeEnd[i] - nodeStart[i]) * WordConfiguration<Word>::WordSize;
		return result;
	}

	void calculateDoublesliceBackwardsRec(size_t i, size_t j, size_t start, Word BA, Word BT, Word BC, Word BG, LengthType sizeLeft, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		assert(sizeLeft > 0);
		bool forceSource = true;
		if (nodeEnd[i]-nodeStart[i] < sizeLeft)
		{
			sizeLeft -= nodeEnd[i]-nodeStart[i];
			for (size_t neighbori = 0; neighbori < inNeighbors[i].size(); neighbori++)
			{
				//todo what happens when two cycles are with 2*w of eachothers?
				//the slice calculation may overwrite a previous correct value with a new wrong value
				//let's assume the cycles are far enough from each others
				//cycling over itself is fine
				assert(inNeighbors[i][neighbori] == start || !notInOrder[inNeighbors[i][neighbori]]);
				auto neighbor = inNeighbors[i][neighbori];
				if (!currentBand[neighbor]) continue;
				calculateDoublesliceBackwardsRec(neighbor, j, start, BA, BT, BC, BG, sizeLeft, currentSlice, previousSlice, currentBand, previousBand);
				forceSource = false;
			}
		}
		calculateNode(i, j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, forceSource);
	}

	void calculateDoublesliceBackwards(size_t j, Word BA, Word BT, Word BC, Word BG, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		if (firstInOrder == 0) return;
		for (size_t i = firstInOrder-1; i > 0; i--)
		{
			std::set<size_t> currentCalculables;
			std::vector<size_t> currentCalculableOrder;
			assert(notInOrder[i]);
			if (!currentBand[i]) continue;
			calculateDoublesliceBackwardsRec(i, j, i, BA, BT, BC, BG, WordConfiguration<Word>::WordSize*2, currentSlice, previousSlice, currentBand, previousBand);
		}
	}

	MatrixSlice getBitvectorSliceScoresAndFinalPosition(const std::string& sequence, int dynamicWidth, std::vector<std::vector<bool>>& startBand, LengthType dynamicRowStart) const
	{
		MatrixSlice result;
		result.cellsProcessed = 0;
		result.finalMinScore = 0;
		result.finalMinScoreColumn = 0;
		result.minScorePerWordSlice.emplace_back(0);

		std::vector<WordSlice> slice1;
		std::vector<WordSlice> slice2;
		slice1.resize(nodeSequences.size());
		slice2.resize(nodeSequences.size());

		std::vector<WordSlice>& currentSlice = slice1;
		std::vector<WordSlice>& previousSlice = slice2;

		std::vector<ScoreType> nodeMinScores;
		nodeMinScores.resize(nodeStart.size(), std::numeric_limits<ScoreType>::max());

		LengthType previousMinimumIndex;
		std::vector<bool> currentBand;
		std::vector<bool> previousBand;
		assert(startBand.size() > 0);
		assert(startBand[0].size() == nodeStart.size());
		currentBand.resize(nodeStart.size());
		previousBand = startBand[0];

		for (size_t j = 0; j < sequence.size(); j += WordConfiguration<Word>::WordSize)
		{
			ScoreType currentMinimumScore = std::numeric_limits<ScoreType>::max();
			LengthType currentMinimumIndex = 0;
			//preprocessed bitvectors for character equality
			Word BA = WordConfiguration<Word>::AllZeros;
			Word BT = WordConfiguration<Word>::AllZeros;
			Word BC = WordConfiguration<Word>::AllZeros;
			Word BG = WordConfiguration<Word>::AllZeros;
			for (int i = 0; i < WordConfiguration<Word>::WordSize && j+i < sequence.size(); i++)
			{
				Word mask = ((Word)1) << i;
				switch(sequence[j+i])
				{
					case 'A':
					case 'a':
					BA |= mask;
					break;
					case 'T':
					case 't':
					BT |= mask;
					break;
					case 'C':
					case 'c':
					BC |= mask;
					break;
					case 'G':
					case 'g':
					BG |= mask;
					break;
					case 'N':
					case 'n':
					BA |= mask;
					BC |= mask;
					BT |= mask;
					BG |= mask;
					break;
				}
			}
			size_t slice = j / WordConfiguration<Word>::WordSize;
			if (startBand.size() > slice)
			{
				if (slice > 0) previousBand = currentBand;
				currentBand = startBand[slice];
			}
			else
			{
				std::swap(currentBand, previousBand);
				currentBand.assign(currentBand.size(), false);
				expandBandFromPrevious(currentBand, previousBand, previousSlice, dynamicWidth, nodeMinScores);
			}
			calculateDoublesliceBackwards(j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand);
			for (size_t i = firstInOrder; i < nodeStart.size(); i++)
			{
				nodeMinScores[i] = std::numeric_limits<ScoreType>::max();
				if (!currentBand[i]) continue;
				auto nodeCalc = calculateNode(i, j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, false);
				nodeMinScores[i] = nodeCalc.minScore;
				if (nodeCalc.minScore < currentMinimumScore)
				{
					currentMinimumScore = nodeCalc.minScore;
					currentMinimumIndex = nodeCalc.minScoreIndex;
				}
				result.cellsProcessed += nodeCalc.cellsProcessed;
			}
			for (size_t i = 0; i < firstInOrder; i++)
			{
				nodeMinScores[i] = std::numeric_limits<ScoreType>::max();
				if (!currentBand[i]) continue;
				auto nodeCalc = calculateNode(i, j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, false);
				nodeMinScores[i] = nodeCalc.minScore;
				if (nodeCalc.minScore < currentMinimumScore)
				{
					currentMinimumScore = nodeCalc.minScore;
					currentMinimumIndex = nodeCalc.minScoreIndex;
				}
				result.cellsProcessed += nodeCalc.cellsProcessed;
			}
			std::swap(currentSlice, previousSlice);
			previousMinimumIndex = currentMinimumIndex;
			result.minScorePerWordSlice.emplace_back(currentMinimumScore);
		}
		result.finalMinScoreColumn = previousMinimumIndex;
		result.finalMinScore = result.minScorePerWordSlice.back();
		return result;
	}

	std::tuple<ScoreType, int, std::vector<MatrixPosition>, size_t> getBacktrace(std::string sequence, int dynamicWidth, LengthType dynamicRowStart, std::vector<std::vector<bool>>& startBand) const
	{
		int padding = (WordConfiguration<Word>::WordSize - (sequence.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		for (int i = 0; i < padding; i++)
		{
			sequence += 'N';
		}
		auto slice = getBitvectorSliceScoresAndFinalPosition(sequence, dynamicWidth, startBand, dynamicRowStart);
		auto backtraceresult = backtrace(std::make_pair(slice.finalMinScoreColumn, sequence.size()), sequence, slice.minScorePerWordSlice);
		for (int i = 0; i < padding; i++)
		{
			assert(backtraceresult.back().second >= sequence.size() - padding);
			backtraceresult.pop_back();
		}
		//if there's a mismatch at the last base, the backtrace might be spending one more jump in the padding
		if (backtraceresult.back().second == sequence.size() - padding && nodeSequences[backtraceresult.back().first] != sequence[backtraceresult.back().second])
		{
			backtraceresult.pop_back();
		}
		assert(backtraceresult.back().second == sequence.size() - padding - 1);
		return std::make_tuple(slice.finalMinScore, 0, backtraceresult, slice.cellsProcessed);
	}

	ScoreType gapPenalty(LengthType length) const
	{
		if (length == 0) return 0;
		return gapStartPenalty + gapContinuePenalty * (length - 1);
	}

	ScoreType matchScore(char graph, char sequence) const
	{
		return graph == sequence ? 1 : -1;
	}

	std::vector<bool> notInOrder;
	std::vector<LengthType> nodeStart;
	std::vector<LengthType> nodeEnd;
	std::vector<LengthType> indexToNode;
	std::map<int, LengthType> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::vector<LengthType>> inNeighbors;
	std::vector<std::vector<LengthType>> outNeighbors;
	std::vector<bool> reverse;
	std::string nodeSequences;
	ScoreType gapStartPenalty;
	ScoreType gapContinuePenalty;
	LengthType dummyNodeStart = 0;
	LengthType dummyNodeEnd = 1;
	bool finalized;
	size_t firstInOrder = 0;
};

#endif