// FILE: Data.cpp (part of namespace rlair_multi_clustering)
// CLASS implemented: Data (see Data.h for documentation)

#include <cstdlib>
#include <cstring>
#include <fstream>                  // provides: ifstream, ofstream
#include <sstream>                  // provides: stringstream
#include <cassert>                  // provides: assert
#include "Data.h"
#include "Indexer.h"

using namespace std;

namespace rlair_multi_clustering
{
  Data::Data() {}

  Data::Data(const string dir, const string data_file, const string labels_file)
    : dir(dir), data_file(data_file), labels_file(labels_file) {}

  void Data::load()
  // Library facilities used: fstream, stringstream
  {
	  // Open data file
    ifstream in;
    string pathname = string(dir + data_file);
	  in.open(pathname.c_str());
	  if(in.fail()) {cerr << "error opening " << pathname << endl; exit(1);}

    // Get modes, dimensions, ways, way-mode map, values
    char line[LENGTH];
    in.getline(line, LENGTH);
    stringstream token(line);
    int modes; token >> modes;
    vector<int> dimension(modes, -1);
    for(int i = 0; i != modes; ++i) token >> dimension[i];
    int ways; token >> ways;
    vector<int> way_mode(ways, -1);
    for(int i = 0; i != ways; ++i) token >> way_mode[i];
    token >> values;

    // Initialize matrix
    vector<int> dimensions(ways);
    for(int i = 0; i != ways; ++i) dimensions[i] = dimension[way_mode[i]];
    matrix = HyperMatrix<T>(dimensions);

    // Get data
    in.getline(line, LENGTH);
    while(!in.eof())
    {
      stringstream token(line);
      vector<int> tuple(ways);
      for(int way = 0; way != ways; ++way) token >> tuple[way];
      T value; token >> value;
      matrix[hyper_index(tuple, dimensions)] = value == 0 ? 0 : 1;
      in.getline(line, LENGTH);
    }

    // Close data file
	  in.close();

    // Open labels file
    in.clear();
    pathname = string(dir + labels_file);
	  in.open(pathname.c_str());
	  if(in.fail()) {cerr << "error opening " << pathname << endl; exit(1);}

    // Read in labels
	  cout << endl; // TEST REMOVE
    int way = 0;
    int count = 0;
    labels = labels_t(ways);
    in.getline(line, LENGTH);
    while(!in.eof())
    {
      cout << strcmp(line, "") << endl;//TEST REMOVE
      if(strcmp(line, ""))
      {
        //for(int way = 0; way != ways; ++way)
          //if(way_mode[way] == mode) labels[way].push_back(string(line));
        labels[way].push_back(string(line));
        ++count;
      }
      else {count = 0; ++way;}
      in.getline(line, LENGTH);
    }

    // check
    for(way = 0; way != ways; ++way) {
      cout << labels[way].size() << " = " << dimension[way] << endl;//TEST REMOVE
      assert(int(labels[way].size()) == dimension[way]);
    }

    // close file
    in.close();
  }

  void Data::operator =(const Data& source)
  // Library facilities used: none
  {
    values = source.values;
    labels = source.labels;
    matrix = source.matrix;
  }

  std::vector<int> Data::dimensions() const
  // Library facilities used: none
  {
    return matrix.dimensions;
  }

  int Data::ways() const
  // Library facilities used: none
  {
    return static_cast<int>(matrix.dimensions.size());
  }
}
