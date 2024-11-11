#ifndef ParallelBB_h
#define ParallelBB_h

#include <thread>
#include <vector>
#include <unistd.h>
#include "BooPHF/BooPHF.h"


template <typename Key>
class ParallelBB
{
	typedef boomphf::SingleHashFunctor<Key> hasher_t;
	typedef boomphf::mphf<Key, hasher_t> boophf_t;
public:
	ParallelBB() {}
	ParallelBB(size_t numThreads, const std::vector<Key>& keys)
	{
		build(numThreads, keys);
	}
	void build(size_t _numThreads, const std::vector<Key>& keys)
	{
		numThreads = _numThreads;
		if (numThreads == 1)
		{
			bucketOffset.resize(2);
			bucketOffset[0] = 0;
			bucketOffset[1] = keys.size();
			buckets.resize(1);
			buckets[0] = new boomphf::mphf<Key,hasher_t>(keys.size(), keys, 1, 2, true, false);
			return;
		}
		std::vector<std::thread> builderThreads;
		buckets.resize(numThreads, nullptr);
		bucketOffset.resize(numThreads+1, 0);
		for (size_t i = 0; i < numThreads; i++)
		{
			builderThreads.emplace_back([this, &keys, i](){
				std::vector<Key> keysHere;
				for (auto key : keys)
				{
					if (hash(key) == i) keysHere.push_back(key);
				}
				bucketOffset[i+1] = keysHere.size();
				buckets[i] = new boomphf::mphf<Key,hasher_t>(keysHere.size(), keysHere, 1, 2, true, false);
			});
		}
		for (size_t i = 0; i < numThreads; i++)
		{
			builderThreads[i].join();
		}
		for (size_t i = 1; i < bucketOffset.size(); i++)
		{
			bucketOffset[i] += bucketOffset[i-1];
		}
	}
	ParallelBB(const ParallelBB& other) = delete;
	ParallelBB& operator=(const ParallelBB& other) = delete;
	ParallelBB(ParallelBB&& other) = default;
	ParallelBB& operator=(ParallelBB&& other) = default;
	~ParallelBB()
	{
		for (auto b : buckets) delete b;
	}
	size_t lookup(const Key& key) const
	{
		size_t bucket = hash(key);
		auto result = buckets[bucket]->lookup(key);
		if (result == ULLONG_MAX) return -1;
		return result + bucketOffset[bucket];
	}
	size_t size() const
	{
		return bucketOffset.back();
	}
private:
	size_t hash(const Key& key) const
	{
		return key % numThreads;
	}
	size_t numThreads;
	std::vector<boophf_t*> buckets;
	std::vector<size_t> bucketOffset;
};

#endif
