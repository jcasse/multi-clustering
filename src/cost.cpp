#include <cmath>                // provides: log
#include "Multiclustering.h"

using namespace std;

namespace rlair_multi_clustering
{
  double Multiclustering::cost()
  // Library facilities used: none
  {return model_encoding_cost() + data_encoding_cost();}

  double Multiclustering::model_encoding_cost()
  // Library facilities used: log
  {
    // (1) send data matrix dimensions (e.g., M and N) using
    // (e.g., log*(m) + log*(n) bits)
    // (generalized to multiple dimenisons)
    // NOTE: this does not change between different clusterings, so it does
    // not factor into the optimization

    // (2) cost of encoding row/column cluster assignments
    // M lg(K) + N lg(L)
    // (generalized to multiple dimenisions)
    int ways = data->ways();
    double assignment_encoding = 0;
    for(int way = 0; way != ways; ++way)
      assignment_encoding += 
      data->matrix.dimensions[way] * log(double(clusterings[way].size()));

    // (3) cost of encoding the type of each block
    // type := distribution of values in block (send number of 1s in blocks)
    // cost per block := (#values - 1) * lg(block_size + 1)
    // NOTE: can make this encoding shorter
    double type_encoding = 0;
    size_t values = data->values - 1;
    Indexer::mask_t mask(ways);
    Indexer::dimensions_t dimensions = blocking_dimensions();
    Indexer indexer(dimensions, mask);
    while(!indexer.end())
    {
      int size = block_size(indexer.get_tuple());
      type_encoding += values * log(double(size) + 1);
      indexer.forward();
    }

    // description complexity cost
    return assignment_encoding + type_encoding;
  }

  double Multiclustering::data_encoding_cost()
  // Library facilities used: log
  {
    double data_encoding = 0;
    Indexer::mask_t mask(data->ways());
    Indexer::dimensions_t dimensions = blocking_dimensions();
    Indexer indexer(dimensions, mask);
    while(!indexer.end())
    {
      data_encoding += block_cost(indexer.get_tuple());
      indexer.forward();
    }
    return data_encoding;
  }

  double Multiclustering::block_cost(const Indexer::tuple_t & tuple) const
  // Library facilities used: none
  {
    counts_t block_counts;
    get_block_counts(block_counts, tuple);
    return hoffman_coding(block_counts, block_frequencies(tuple));
  }

  frequencies_t Multiclustering::block_frequencies
  (const Indexer::tuple_t & tuple) const
  // Library facilities used: none
  {
    int total = block_size(tuple);
    counts_t value_counts;
    get_block_counts(value_counts, tuple);
    frequencies_t value_frequencies(data->values);
    for(int i = 0; i != data->values; ++i)
      value_frequencies[i] = frequency(value_counts[i], total);
    return value_frequencies;
  }

  void Multiclustering::get_block_counts
  (counts_t & counts, const Indexer::tuple_t & tuple) const
  // Library facilities used: none
  {
    counts = counts_t(data->values);
    Indexer::indexes_t indexes(data->ways());
    for(int way = 0; way != data->ways(); ++way)
      indexes[way] = clusterings[way][tuple[way]];
    Indexer::mask_t mask(data->ways());
    Indexer indexer(indexes, mask);
    while(!indexer.end())
    {
      vector<int> dimensions = data->dimensions();
      ++counts[data->matrix[hyper_index(indexer.get_tuple(), dimensions)]];
      indexer.forward();
    }
  }

  double hoffman_coding(const vector<int> & counts)
  {
    int total = 0;
    for(size_t value = 0; value != counts.size(); ++value)
      total += counts[value];
    double cost = 0;
    for(size_t v = 0; v != counts.size(); ++v)
      cost += hoffman_coding(counts[v], frequency(counts[v], total));
    return cost;
  }

  double hoffman_coding
  (const vector<int> & counts, const vector<double> & frequencies)
  // Library facilities used: assert
  {
    assert(counts.size() == frequencies.size());
    double cost = 0;
    for(size_t i = 0; i != counts.size(); ++i)
      cost += hoffman_coding(counts[i], frequencies[i]);
    return cost;
  }
 
  double hoffman_coding
  (const vector<int> & unit_counts, const vector<int> & block_counts)
  // Library facilities used: assert
  {
    assert(unit_counts.size() == block_counts.size());
    int total = 0;
    for(size_t i = 0; i != block_counts.size(); ++i) total += block_counts[i];
    double cost = 0;
    for(size_t i = 0; i != unit_counts.size(); ++i)
      cost +=
        hoffman_coding(unit_counts[i], frequency(block_counts[i], total));
    return cost;
  }
 
  double hoffman_coding(int count, double frequency)
  // Library facilities used: assert, log
  {
    if(count == 0) return 0;
    if(frequency == 0) return DBL_MAX;
    return count * -log(frequency);
  }

  double frequency(int count, int total)
  // Library facilities used: none
  {return count / double(total);}
}
