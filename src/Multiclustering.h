// FILE: Multiclustering.h
// CLASS PROVIDED: Multiclustering (part of the namespace rlair_multi_clustering)

#ifndef RLAIR_MULTI_CLUSTERING_MULTICLUSTERING
#define RLAIR_MULTI_CLUSTERING_MULTICLUSTERING

// FILES
#include <iostream>                 // provides: ostream
#include <iomanip>                  // provides: setw, fixed, setprecision
#include <climits>                  // provides: INT_MAX
#include <cfloat>                   // provides: DBL_MAX

#include "Data.h"
#include "Indexer.h"

#define BOOST_FILESYSTEM_VERSION 3

//#include "boost/shared_ptr.hpp"     // provides: shared_ptr

namespace rlair_multi_clustering
{
  typedef std::vector<int> counts_t;
  typedef std::vector<int> cluster_t;
  typedef std::vector<double> frequencies_t;
  typedef std::vector<cluster_t> clustering_t;
  typedef std::vector<cluster_t> multicluster_t;

  class Multiclustering
  {
  public:
    // CONSTRUCTORS and DESTRUCTOR
    Multiclustering();
    Multiclustering(Data * data, std::ostream * lout);
    Multiclustering(Data * data, std::ostream * lout, std::vector<int> & clusters);
    Multiclustering(const Multiclustering & source);

    // MODIFICATION MEMBER FUNCTIONS
    void copy(const Multiclustering & source);
    // Precondition: none
    // Postcondition: *this is a copy of source
    void initialize();
    // Precondition: none
    // Postcondition: *this has valid initial data
    void initialize(std::vector<int> & clusters);
    // Precondition: none
    // Postcondition: *this has valid initial data and clusterings of clusters
    Multiclustering & operator=(const Multiclustering & source);
    // Precondition: none
    // Postcondition: *this == source
    void make_multiclustering(std::vector<int> & clusters);
    // Precondition: nonde
    // Postcondition: clusterings have clusters per mode
    bool optimize(int way);
    // Precondition: way < matrix ways
    // Postcondition: all units in way have been placed in best clusters
    bool add_cluster(int way);
    // Precondition: way is a valid data matrix way
    // Postcondition: multiclustering has one additional cluster in way

    // CONSTANT MEMBER FUNCTIONS
    double cost();
    // Precondition: data, clusterings, and blocks have been initialized
    // Postcondition: Return value is data encoding cost
    double model_encoding_cost();
    // Precondition: data, clusterings, and blocks have been initialized
    // Postcondition: Return value is model description length
    double data_encoding_cost();
    // Precondition: blocks has been initialized
    // Postcondition: Return value is data description length
    void print_2D_slice
    (const std::vector<int> & dimension, const std::string & file) const;
    void print_2D_slice
    (const std::vector<int> & dimension, std::ostream & out) const;
    void print_model_2D
    (const std::vector<int> & dimension, const std::string & file) const;
    void print_model_2D
    (const std::vector<int> & dimension, std::ostream & out) const;
    void print_blocked_matrix_2D(const std::string & file) const;
    void print_blocked_matrix_2D(std::ostream & out) const;
    void print_clusterings(const std::string & dir) const;
    void print_block_densities(const std::string & dir) const;

    // MEMBER VARIABLES
    Data * data;                                // pointer to data
    std::vector<clustering_t> clusterings;      // clusters of units
    std::ostream * lout;                        // pointer to log file

  private:

    // UTILITY MEMBER FUNCTIONS
    void print_clustered_row
    (int ROW, int cluster, int row, int COL, int cluster_z, int index_z, std::ostream & out, bool model) const;
    // Precondition: cluster is row cluster, row is row in cluster
    // cluster == -1 prints clusteirng heading
    // row == -1 prints top line, row == # rows prints bottom line
    // Postcondition: clustered row belonging to a slice is printed
    int blocking_size() const;
    // Precondition: none
    // Postcondition: Return value is number of blocks
    int blocking_size(int way) const;
    // Precondition: none
    // Postcondition: Return number of blocks in hyper-plane with way fixed
    Indexer blocking_indexer(int way, int cluster);
    // Precondition: way within data ways, cluster within way clustering
    // Postcondition: Return value is indexer for blocks around way cluster
    Indexer block_indexer(multicluster_t & block, int way, int unit_index);
    // Precondition: way < block ways, unit < way cluster size
    // Postcondition: Return value is indexer for block units around way unit
    multicluster_t get_block(const Indexer::tuple_t & tuple);
    // Precondition: tuple is valid
    // Postcondition: Return value is multicluster indexed by tuple

    void get_unit_signature
    (std::vector<counts_t> & signature, int way, int cluster, int unit_index);
    // Precondition: way < ways, cluster < way clusters, unit < cluster units
    // Postcondition: Return value is unit's value counts for each block
    bool optimize(int way, int unit);
    // Precondition: way < matrix ways, unit < way units
    // Postcondition: way unit is placed in optimum way cluster
    // faster (04/22/12)
    bool optimize(int way, int old_cluster, int index,
      std::vector<std::vector<counts_t> > & units_signatures,
      std::vector<std::vector<counts_t> > & clusters_signatures,
      std::vector<int> & new_assignments);
    // Precondition: way < matrix ways, unit < way units
    // Postcondition: way unit is placed in optimum way cluster
    int split_cluster(int way);
    // Precondition: way is a valid data matrix way
    // Postcondition: Return value is index of cluster to split
    double cluster_cost(const std::vector<counts_t> & block_counts);
    // Precondition: block_counts for all blocks
    // Postcondition: Return value is cost of blocks
    double cluster_cost(int way, int cluster);
    // Precondition: way is a valid data matrix way, cluster < way clusters
    // Postcondition: Return value is cluster cost

    Indexer::dimensions_t blocking_dimensions() const;
    // Precondition: none
    // Postcondition: Return value is dimensions of the blocking (e.g., K x L)
    int block_size(const Indexer::tuple_t & tuple) const;
    // Precondition: none
    // Postcondition: Return value is number of units in block
    double block_cost(const Indexer::tuple_t & tuple) const;
    // Precondition: none
    // Postcondition: Return value is encoding cost of block
    frequencies_t block_frequencies(const Indexer::tuple_t & tuple) const;
    // Precondition: none
    // Postcondition: Return value is frequency of each value in block
    void get_block_counts
    (counts_t & counts, const Indexer::tuple_t & tuple) const;
    // Precondition: none
    // Postcondition: Return value is number of units in block
  };

  double hoffman_coding
  (const std::vector<int> & counts);
  // Precondition: counts is the number of times each value appears
  // Postcondition: Return value is Hoffman encoding cost in nats
  double hoffman_coding
  (const std::vector<int> & counts, const std::vector<double> & frequencies);
  // Precondition: counts is the number of times a value appears, frequencies
  // is the frequency of the value appearing relative to block size
  // Postcondition: Return value is Hoffman encoding cost in nats
  double hoffman_coding
  (const std::vector<int> & unit_counts, const std::vector<int> & block_counts);
  // Precondition: counts is the number of times a value appears, frequencies
  // is the frequency of the value appearing relative to block size
  // Postcondition: Return value is Hoffman encoding cost in nats
  double hoffman_coding(int count, double frequency);
  // Precondition: count is the number of occurences of an item, and frequency
  // its frequency relative to the total number of items in a block
  // Postcondition: Return value is the hoffman conding cost in nats
  double frequency(int count, int total);
  // Precondition: count is the number of occurrences of a type of item, and
  // total is the total number of all types of items
  // Postcondition: Return value is the frequency (> 0) of the item
}

#endif
