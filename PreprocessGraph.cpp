#include <fstream>
#include <tuple>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "vg.pb.h"
#include "CommonUtils.h"
#include "BigraphToDigraph.h"

//sometimes the graph is really large and doing DFS recursively runs into a stack overflow.
//use our own stack in the heap
//see https://stackoverflow.com/questions/159590/way-to-go-from-recursion-to-iteration
class DFSStack
{
public:
	DFSStack(int i, int state) : i(i), state(state) {};
	int i;
	int state;
};

void topological_sort_using_DFS_stackless(const std::vector<std::vector<size_t>>& graph, std::vector<bool>& explored, std::set<size_t>& mfvs, std::vector<size_t>& currentStack, size_t i, std::vector<size_t>& sorted, size_t& t)
{
	std::vector<DFSStack> stack;
	stack.emplace_back(i, 1);
	while (!stack.empty())
	{
		auto top = stack.back();
		stack.pop_back();
		i = top.i;
		switch(top.state)
		{
			case 1:
				for (size_t k = 0; k < currentStack.size(); k++)
				{
					if (currentStack[k] == i)
					{
						mfvs.insert(i);
						explored[i] = true;
						continue;
					}
				}
				if (explored[i]) continue;
				currentStack.push_back(i);
				explored[i] = true;
				stack.emplace_back(i, 2);
				for (size_t j = 0; j < graph[i].size(); ++j)
				{
					stack.emplace_back(graph[i][j], 1);
				}
				break;

			case 2:
				assert(currentStack.size() > 0);
				assert(currentStack.back() == i);
				currentStack.pop_back();
				--t;
				sorted[t] = i;
				break;

			default:
				assert(false);
		}
	}
}

void topological_sort_using_DFS_loop(const std::vector<std::vector<size_t>>& graph, std::vector<size_t>& sorted, std::set<size_t>& mfvs)
{
	assert(graph.size() == sorted.size());
	std::vector<bool> explored(graph.size(), false);
	size_t t = graph.size();

	for (size_t i = 0; i < graph.size(); ++i)
	{
		if (explored[i] == false)
		{
			std::vector<size_t> currentStack;
			topological_sort_using_DFS_stackless(graph, explored, mfvs, currentStack, i, sorted, t);
		}
	}

	assert(t == 0);
	std::vector<size_t> newIndex;
	std::set<size_t> resultSet;
	newIndex.resize(sorted.size(), 1);
	newIndex[0] = 0;
	for (auto x : mfvs)
	{
		newIndex[x]--;
	}
	for (size_t i = 1; i < sorted.size(); i++)
	{
		newIndex[i] = newIndex[i] + newIndex[i-1];
	}
	if (mfvs.size() > 0)
	{
		int removed = 0;
		for (size_t i = sorted.size()-1; i < sorted.size(); i--)
		{
			if (mfvs.count(sorted[i]) > 0)
			{
				sorted.erase(sorted.begin()+i);
				removed++;
			}
			else
			{
				sorted[i] = newIndex[sorted[i]];
				assert(sorted[i] >= 0);
				assert(sorted[i] < graph.size() - mfvs.size());
				assert(resultSet.count(sorted[i]) == 0);
				resultSet.insert(sorted[i]);
			}
		}
		assert(removed == mfvs.size());
	}
	assert(resultSet.size() == sorted.size());

	return;
}

std::pair<std::vector<size_t>, std::vector<size_t>> topologicalSort(const DirectedGraph& digraph)
{
	std::vector<std::vector<size_t>> graph;
	graph.resize(digraph.nodes.size());

	for (size_t i = 0; i < digraph.edges.size(); i++)
	{
		assert(digraph.edges[i].fromIndex < graph.size());
		assert(digraph.edges[i].toIndex < graph.size());
		graph[digraph.edges[i].fromIndex].push_back(digraph.edges[i].toIndex);
	}

	std::vector<size_t> sorted(digraph.nodes.size(), 0);
	std::set<size_t> mfvs;
	topological_sort_using_DFS_loop(graph, sorted, mfvs);
	std::vector<size_t> mfvsVec;
	mfvsVec.insert(mfvsVec.end(), mfvs.begin(), mfvs.end());
	return std::make_pair(mfvsVec, sorted);
}

int main(int argc, char** argv)
{
	auto graph = CommonUtils::LoadVGGraph(argv[1]);
	DirectedGraph digraph {graph};
	auto result = topologicalSort(digraph);
	std::ofstream mfvs {argv[2]};
	boost::archive::text_oarchive mfvsoa(mfvs);
	mfvsoa << result.first;
	std::ofstream order {argv[3]};
	boost::archive::text_oarchive orderoa(order);
	orderoa << result.second;
}