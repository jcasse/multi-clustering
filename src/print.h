#include <fstream>              // provides: ifstream, ofstream
#include "Indexer.h"
#include "Multiclustering.h"

using namespace std;

namespace rlair_multi_clustering
{
  void Multiclustering::print_2D_slice
  (const vector<int> & dimension, const string & file) const
  {
    // open file
    ofstream out;
    out.open(file.c_str());
    if(out.fail()) {cout << "error opening " << file << endl; exit(1);}
    out << fixed << setprecision(2);

    // print co-clustered matrix
    print_2D_slice(dimension, out);

    // close file
    out.close();
  }

  void Multiclustering::print_2D_slice
  (const vector<int> & dimension, ostream & out) const
  {
    int ways = int(data->ways());

    // get plane to print
    vector<int> plane;
    for(int way = 0; way != ways; ++way)
      if(dimension[way] == -1) plane.push_back(way);

    assert(int(plane.size()) == 2);

    const int ROW = plane[0];
    const int COL = plane[1];

    assert(ROW < COL);

    const int MAX_DIMENSION = 25;
    if(data->matrix.dimensions[ROW] > MAX_DIMENSION)
    {cout << "number of rows > " << MAX_DIMENSION << endl; return;}
    if(data->matrix.dimensions[COL] > MAX_DIMENSION)
    {cout << "number of columns > " << MAX_DIMENSION << endl; return;}

    // get cluster_z and index_z
    int way_z = -1;
    int index_z = -1;
    int cluster_z = -1;
    if(ways == 3)
    {
    if(ROW + COL == 1) way_z = 2;
    if(ROW + COL == 2) way_z = 1;
    if(ROW + COL == 3) way_z = 0;
    for(int cluster = 0; cluster != int(clusterings[way_z].size()); ++cluster)
    {
      int indexes = int(clusterings[way_z][cluster].size());
      for(int index = 0; index != indexes; ++index)
      {
        if(clusterings[way_z][cluster][index] == dimension[way_z])
        {
          index_z = index;
          cluster_z = cluster;
        }
      }
    }
    }

    // initialize variables
    int rows = 1;
    const int HEADING = -1;
    const int heading = 0;
    const int topline = -1;
    const int & botline = rows;

    // print heading
    print_clustered_row
    (ROW, HEADING, topline, COL, cluster_z, index_z, out, false);
    print_clustered_row
    (ROW, HEADING, heading, COL, cluster_z, index_z, out, false);
    print_clustered_row
    (ROW, HEADING, botline, COL, cluster_z, index_z, out, false);

    // print row clusters
    int clusters = int(clusterings[ROW].size());
    for(int cluster = 0; cluster != clusters; ++cluster)
    {
      // print top line
      print_clustered_row
      (ROW, cluster, topline, COL, cluster_z, index_z, out, false);

      // print rows
      rows = int(clusterings[ROW][cluster].size());
      for(int row = 0; row != rows; ++row)
      {
        print_clustered_row
        (ROW, cluster, row, COL, cluster_z, index_z, out, false);
      }

      // print bottom line
      print_clustered_row
      (ROW, cluster, botline, COL, cluster_z, index_z, out, false);
    }
  }

  void Multiclustering::print_model_2D
  (const vector<int> & dimension, const string & file) const
  {
    // open file
    ofstream out;
    out.open(file.c_str());
    if(out.fail()) {cout << "error opening " << file << endl; exit(1);}
    out << fixed << setprecision(2);

    // print co-clustered matrix
    print_model_2D(dimension, out);

    // close file
    out.close();
  }
  
  void Multiclustering::print_model_2D
  (const vector<int> & dimension, ostream & out) const
  {
    int ways = int(data->ways());

    // get plane to print
    vector<int> plane;
    for(int way = 0; way != ways; ++way)
      if(dimension[way] == -1) plane.push_back(way);

    assert(int(plane.size()) == 2);

    const int ROW = plane[0];
    const int COL = plane[1];

    assert(ROW < COL);

    // get cluster_z and index_z
    int way_z = -1;
    int index_z = -1;
    int cluster_z = -1;
    if(ways == 3)
    {
    if(ROW + COL == 1) way_z = 2;
    if(ROW + COL == 2) way_z = 1;
    if(ROW + COL == 3) way_z = 0;
    for(int cluster = 0; cluster != int(clusterings[way_z].size()); ++cluster)
    {
      int indexes = int(clusterings[way_z][cluster].size());
      for(int index = 0; index != indexes; ++index)
      {
        if(clusterings[way_z][cluster][index] == dimension[way_z])
        {
          index_z = index;
          cluster_z = cluster;
        }
      }
    }
    }

    // initialize variables
    int rows = 1;
    const int HEADING = -1;
    const int heading = 0;
    const int topline = -1;
    const int & botline = rows;
    const int model = -2;

    // print heading
    print_clustered_row
    (ROW, HEADING, topline, COL, cluster_z, index_z, out, true);
    print_clustered_row
    (ROW, HEADING, heading, COL, cluster_z, index_z, out, true);
    print_clustered_row
    (ROW, HEADING, botline, COL, cluster_z, index_z, out, true);

    // print row clusters
    int clusters = int(clusterings[ROW].size());
    for(int cluster = 0; cluster != clusters; ++cluster)
    {
      // print top line
      print_clustered_row
      (ROW, cluster, topline, COL, cluster_z, index_z, out, true);

      // print rows
      print_clustered_row
      (ROW, cluster, model, COL, cluster_z, index_z, out, true);

      // print bottom line
      print_clustered_row
      (ROW, cluster, botline, COL, cluster_z, index_z, out, true);
    }
  }

  void Multiclustering::print_clustered_row(int ROW, int cluster, int row,
  int COL, int cluster_z, int index_z, ostream & out, bool model)
  const
  // Library facilities used: none
  {
    // graphics constants
    const char ltangle = char(218);
    const char rtangle = char(191);
    const char lbangle = char(192);
    const char rbangle = char(217);
    const char vertbar = char(179);
    const char horzbar = char(196);
    const char hzspace = char(32);

    // 2D constants
    const int HEADING = -1;

    // initialize variables
    int rows = 1;
    if(cluster != HEADING) rows = int(clusterings[ROW][cluster].size());
    if(model) rows = 1;
    const int topline = -1;
    const int & botline = rows;
    int row_unit = -1;
    if(cluster != HEADING && topline < row && row < botline)
      row_unit = clusterings[ROW][cluster][row];
    int way_z = -1;
    if(ROW + COL == 1) way_z = 2;
    if(ROW + COL == 2) way_z = 1;
    if(ROW + COL == 3) way_z = 0;

    // print left heading
    if(cluster == HEADING)
    {
      if(cluster_z < 0 || row == topline || row == botline)
      out << hzspace << hzspace << hzspace << hzspace;
      else if(index_z != -1)
        out << setw(3) << clusterings[way_z][cluster_z][index_z] << hzspace;
    }
    else if(row == topline) out << ltangle << horzbar << horzbar << rtangle;
    else if(row == botline) out << lbangle << horzbar << horzbar << rbangle;
    else
    {
      if(model) out << vertbar << setw(2) << char(cluster + 65) << vertbar;
      else out << vertbar << setw(2) << row_unit << vertbar;
    }
    out << hzspace;

    // print clustered row
    int clusters = static_cast<int>(clusterings[COL].size());
    for(int column_cluster = 0; column_cluster != clusters; ++column_cluster)
    {
      // print left end symbol
      if(row == topline) out << ltangle;
      else if(row == botline) out << lbangle;
      else out << vertbar;

      // print contents
      if(model)
      {
        if(row == topline) out << horzbar << horzbar << horzbar << horzbar;
        else if(row == botline) out << horzbar << horzbar << horzbar << horzbar;
        else if(cluster == HEADING)
          out << "  " << char(column_cluster + 65) << " ";
        else
        {
          // compute block value frequency
          char frequency[6];
          Indexer::tuple_t tuple;
          if(cluster_z < 0)
          {
            tuple = Indexer::tuple_t(2, cluster);
            tuple[COL] = column_cluster;
          }
          else
          {
            tuple = Indexer::tuple_t(3, cluster);
            tuple[COL] = column_cluster;
            if(ROW + COL == 1) tuple[2] = cluster_z;
            if(ROW + COL == 2) tuple[1] = cluster_z;
            if(ROW + COL == 3) tuple[0] = cluster_z;
          }
          frequencies_t frequencies = block_frequencies(tuple);
          sprintf(frequency, "%0.3f", frequencies[1]);
          out << frequency;
        }
      }
      else
      {
        int units = static_cast<int>(clusterings[COL][column_cluster].size());
        for(int unit_index = 0; unit_index != units; ++unit_index)
        {
          if(row == topline) out << horzbar << horzbar;
          else if(row == botline) out << horzbar << horzbar;
          else if(cluster == HEADING)
            out << setw(2) << clusterings[COL][column_cluster][unit_index];
          else
          {
            Indexer::tuple_t tuple;
            if(cluster_z < 0)
            {
              tuple = Indexer::tuple_t(2, row_unit);
              tuple[COL] = clusterings[COL][column_cluster][unit_index];
            }
            else
            {
              tuple = Indexer::tuple_t(3, row_unit);
              tuple[COL] = clusterings[COL][column_cluster][unit_index];
              if(ROW + COL == 1) tuple[2] = clusterings[2][cluster_z][index_z];
              if(ROW + COL == 2) tuple[1] = clusterings[1][cluster_z][index_z];
              if(ROW + COL == 3) tuple[0] = clusterings[0][cluster_z][index_z];
            }
            out << setw(2) <<
            data->matrix.data[hyper_index(tuple, data->matrix.dimensions)];
          }
        }
      }

      // print final space
      if(row == topline) out << horzbar;
      else if(row == botline) out << horzbar;
      else if(cluster == HEADING || !model) out << hzspace;

      // right end symbol
      if(row == topline) out << rtangle;
      else if(row == botline) out << rbangle;
      else out << vertbar;
    }
    out << endl;
  }

  void Multiclustering::print_clusterings(const string & dir) const
  // Library facilities used: none
  {
    int ways = data->ways();
    for(int way = 0; way < ways; way++)
    {
      // open file
      char id[LENGTH];
      sprintf(id, "%d", way);
      string file = string(dir + "clustering_" + id + ".txt");
      ofstream out;
      out.open(file.c_str());
      if(out.fail()) {cout << "error opening " << file << endl; exit(1);}

      // print clustering
      int clusters = int(clusterings[way].size());
      for(int cluster = 0; cluster < clusters; cluster++)
      {
        out << "cluster " << cluster << endl;
        int units = int(clusterings[way][cluster].size());
        for(int unit = 0; unit < units; unit++)
        {
          out << setw(6) << clusterings[way][cluster][unit] << " ";
          out << data->labels[way][clusterings[way][cluster][unit]] << endl;
        }
        if(cluster < clusters - 1) out << endl;
      }

      // close file
      out.close();
    }
  }

  void Multiclustering::print_block_densities(const string & file) const
  // Library facilities used: none
  {
    // open file
    ofstream out;
    out.open(file.c_str());
    if(out.fail()) {cout << "error opening " << file << endl; exit(1);}
    out << fixed << setprecision(2);

    // create indexes
    int ways = data->ways();
    Indexer::indexes_t indexes(ways);
    for(int i = 0; i != ways; ++i)
      for(int j = 0; j != int(clusterings[i].size()); ++j)
        indexes[i].push_back(j);

    // set tuple
    vector<int> tuple(ways);

    // set mask for way
    Indexer::mask_t mask(ways);

    // create indexer
    Indexer indexer(indexes, tuple, mask);

    while(!indexer.end())
    {
      Indexer::tuple_t tuple = indexer.get_tuple();
      frequencies_t freqs = block_frequencies(tuple);
      for(size_t i = 0; i != tuple.size(); ++i) out << tuple[i] << " ";
      for(size_t i = 0; i != freqs.size(); ++i) out << freqs[i] << " ";
      out << endl;
      indexer.forward();
    }

    // close file
    out.close();
  }

  void Multiclustering::print_blocked_matrix_2D(const string & file) const
  // Library facilities used: none
  {
    // open file
    ofstream out;
    out.open(file.c_str());
    if(out.fail()) {cout << "error opening " << file << endl; exit(1);}
    out << fixed << setprecision(2);

    // print co-clustered matrix
    print_blocked_matrix_2D(out);

    // close file
    out.close();
  }

  void Multiclustering::print_blocked_matrix_2D(ostream& out) const
  // Library facilities used: none
  {
    int ways = data->ways();
    vector<int> tuple(ways, 0);
    int row_clusters = int(clusterings[0].size());
    for(int row_cluster = 0; row_cluster < row_clusters; row_cluster++) {
      int row_cluster_elements = int(clusterings[0][row_cluster].size());
      for(int row_cluster_element = 0; row_cluster_element < row_cluster_elements; row_cluster_element++) {
        tuple[0] = clusterings[0][row_cluster][row_cluster_element];
        int column_clusters = int(clusterings[1].size());
        for(int column_cluster = 0; column_cluster < column_clusters; column_cluster++) {
          int column_cluster_elements = int(clusterings[1][column_cluster].size());
          for(int column_cluster_element = 0; column_cluster_element < column_cluster_elements; column_cluster_element++) {
            tuple[1] = clusterings[1][column_cluster][column_cluster_element];
            if(data->matrix[hyper_index(tuple, data->matrix.dimensions)] == 0)
              out << 0;
            else out << 1;
          }
          out << char(32);
        }
        out << endl;
      }
      if(row_cluster < row_clusters - 1) out << endl;
    }
  }
}
