// FILE: Data.h
// CLASS PROVIDED: Data (part of the namespace rlair_multi_clustering)

#ifndef RLAIR_MULTI_CLUSTERING_DATA
#define RLAIR_MULTI_CLUSTERING_DATA

#define LENGTH 256

#include <iostream>                 // provides: ostream
#include <string>                   // provides: string
#include "HyperMatrix.h"

namespace rlair_multi_clustering
{
  class Data
  {
  public:
    typedef int T;
    typedef std::vector<std::vector<std::string> > labels_t;

    // CONSTRUCTORS and DESTRUCTOR
    Data();
    Data(const std::string dir, const std::string data_file, const std::string labels_file);

    // MODIFICATION MEMBER FUNCTIONS
    void load();
    // Precondition: none
    // Postcondition: matrix and labels have data from files
    void operator =(const Data & source);
    // Precondition: none
    // Postcondition: *this == source

    // CONSTANT MEMBER FUNCTIONS
    std::vector<int> dimensions() const;
    // Precondition: none
    // Postcondition: Return value is vector of way lengths
    int ways() const;
    // Precondition: none
    // Postcondition: Return value is number of matrix ways

    // MEMBER VARIABLES
    int values;                       // number distinct values
    labels_t labels;                  // labels for each way unit
    HyperMatrix<T> matrix;            // data matrix

  private:
    const std::string dir;
    const std::string data_file;
    const std::string labels_file;
  };
}

#endif
