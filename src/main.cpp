// FILE main.cpp
// Driver program for the rlair_multi_clustering package

// FILES
#include <cstdlib>
#include <cstring>
#include <fstream>                  // provides: ifstream, ofstream
#include <sstream>                  // provides: stringstream
#include <ctime>                    // provides: time, clock, timeinfo, mktime, CLOCKS_PER_SEC
#include <climits>                  // provides: INT_MAX
#include <iomanip>                  // provides: setw, fixed, setprecision
#include <cerrno>                   // provides: errno

// FILES (directory creation in WIN32 / LINUX)
#ifdef WIN32
#include "boost/filesystem.hpp"     // provides: create_directory
#endif
#ifndef WIN32
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "Data.h"
#include "Multiclustering.h"
#include "print.h"

using namespace std;
using namespace rlair_multi_clustering;

// CONSTANTS
const string LOG_FILE("log.txt");
const string DATA_FILE("data.txt");
const string LABELS_FILE("labels.txt");
const string MATRIX_FILE("matrix.txt");
const string BLOCKED_MATRIX_FILE("blocked_matrix.txt");
const string BLOCK_MODEL_FILE("block_model.txt");
const string BLOCK_DENSITIES_FILE("densities.txt");

string bytes(size_t size)
// Library facilities used: none
{
  stringstream ss;
  if(size < 1024) ss << size << " bytes";
  else if(size < 1048576) ss << size / 1024 << " KB";
  else if(size < 1073741824) ss << size / 1048576 << " MB";
  else ss << size << " GB";
  return ss.str();
}

void load_data(Data & data)
// Library facilities used: none
{
  cout << "loading data . . . ";
  time_t start = time(NULL);
  data.load();
  time_t finish = time(NULL);
  cout << bytes(data.matrix.size());
  cout << "done " << finish - start << " seconds" << endl << endl;
}

string timestamp()
// Library facilities used: time, timeinfo, mktime, sprintf, strcpy, strcat
{
  // Get time info
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  mktime(timeinfo);

  // Construct timestamp: YEAR MONTH DAY HOUR MINUTE SECOND
  // example: May 30, 2011 at 5:03:41pm = 20110530170341
  stringstream ss;
  ss << setfill ('0') << setw (4) << timeinfo->tm_year + 1900;
  ss << setfill ('0') << setw (2) << timeinfo->tm_mon + 1;
  ss << setfill ('0') << setw (2) << timeinfo->tm_mday;
  ss << setfill ('0') << setw (2) << timeinfo->tm_hour;
  ss << setfill ('0') << setw (2) << timeinfo->tm_min;
  ss << setfill ('0') << setw (2) << timeinfo->tm_sec;

  return ss.str();
}

string create_output_dir(const string & input_dir, const string & prog_name)
// Library facilities used: boost::filesystem::create_directory
{
  // Set output directory
  string output_dir
    (input_dir + prog_name + "/" + prog_name + "_" + timestamp() + "/");

  cout << "output dir " << output_dir << endl << endl;

  // Create output directory
#ifdef WIN32
  if(!boost::filesystem::create_directory(output_dir.c_str()))
  {
    cout << "error: could not create directory " << output_dir << endl;
    exit(1);
  }
#endif
#ifndef WIN32
  if(mkdir(output_dir.c_str(), 0760) == -1)//creating a directory
  {
    cerr << "Error :  " << strerror(errno) << endl;
    exit(1);
  }
#endif

  return output_dir;
}

//void symmetric_encode_UN_world_trade_data(Data & data)
//// Library facilities used: none
//{
//  // Encode UN world trade data to make it symmetric
//  // symetric (3), export (1), import (2), zero (0)
//  data.values = 4;
//  int units = data.matrix.dimensions[0];
//  for(int i = 0; i != units; ++i)
//  {
//    vector<int> index_upper(data.ways(), i);
//    vector<int> index_lower(data.ways(), i);
//    for(int j = i; j != units; ++j)
//    {
//      index_upper[1] = j;
//      Data::T upper = data.matrix[hyper_index(index_upper, data.dimensions())];
//      index_lower[0] = j;
//      Data::T lower = data.matrix[hyper_index(index_lower, data.dimensions())];
//      if(lower == 1 && upper == 0)
//        data.matrix[hyper_index(index_upper, data.dimensions())] = 2;
//      if(lower == 0 && upper == 1)
//        data.matrix[hyper_index(index_lower, data.dimensions())] = 2;
//      if(lower == 1 && upper == 1)
//      {
//        data.matrix[hyper_index(index_upper, data.dimensions())] = 3;
//        data.matrix[hyper_index(index_lower, data.dimensions())] = 3;
//      }
//    }
//  }
//}

vector<Data::T> symmetric_encode_data(Data & data)
// Library facilities used: none
{
  vector<Data::T> symmetric_mask(data.matrix.size());
  int units = data.matrix.dimensions[0];
  for(int i = 0; i != units; ++i)
  {
    vector<int> index_upper(data.ways(), i);
    vector<int> index_lower(data.ways(), i);
    for(int j = i; j != units; ++j)
    {
      index_upper[1] = j;
      Data::T upper = data.matrix[hyper_index(index_upper, data.dimensions())];
      index_lower[0] = j;
      Data::T lower = data.matrix[hyper_index(index_lower, data.dimensions())];
      if(lower == 1 && upper == 0)
      {
        int index = hyper_index(index_upper, data.dimensions());
        data.matrix[index] = 1;
        symmetric_mask[index] = 1;
      }
      if(lower == 0 && upper == 1)
      {
        int index = hyper_index(index_lower, data.dimensions());
        data.matrix[index] = 1;
        symmetric_mask[index] = 1;
      }
    }
  }

  return symmetric_mask;
}

double regroup(Multiclustering & local, ostream & lout)
// Library facilities used: none
{
  double new_cost = DBL_MAX;
  double old_cost = local.cost();

  cerr << "\t\t\told cost = " << old_cost << endl;
  lout << "\t\t\told cost = " << old_cost << endl;

  while(new_cost != old_cost)
  {
    old_cost = new_cost;

    for(int way = 0; way != local.data->ways(); ++way)
    {
      cerr << "\t\t\toptimize way " << way << " . . ." << endl;
      lout << "\t\t\toptimize way " << way << " . . ." << endl;

      time_t start_01 = time(NULL);
      local.optimize(way);
      time_t finish_01 = time(NULL);

      cerr << "\t\t\ttime = " << finish_01 - start_01 << " seconds" << endl;
      lout << "\t\t\ttime = " << finish_01 - start_01 << " seconds" << endl;
    }

    new_cost = local.cost();

    cerr << "\t\t\tnew cost = " << new_cost << endl;
    lout << "\t\t\tnew cost = " << new_cost << endl;
  }

  return new_cost;
}

Multiclustering crossassociation_search(Data & data, ostream & lout)
// Library facilities used: none
{
  // initialize multiclustering to 1 cluster per way
  Multiclustering global = Multiclustering(&data, &lout);

  double old_cost = DBL_MAX;
  double new_cost = global.cost();

  cerr << "\told cost = " << new_cost << endl;
  lout << "\told cost = " << new_cost << endl;

  while(new_cost != old_cost)
  {
    old_cost = new_cost;

    // initialize best clustering
    Multiclustering best = global;
    double best_cost = global.cost();

    time_t start_01 = time(NULL);

    // pick best way to increment number of clusters
    for(int way = 0; way != data.ways(); ++way)
    {
      Multiclustering local = global;

      cerr << "\t\tadding cluster in way " << way << " . . ." << endl;
      lout << "\t\tadding cluster in way " << way << " . . ." << endl;

      time_t start_02 = time(NULL);
      local.add_cluster(way);
      time_t finish_02 = time(NULL);

      cerr << "\t\ttime = " << finish_02 - start_02 << " seconds" << endl;
      lout << "\t\ttime = " << finish_02 - start_02 << " seconds" << endl;

      cerr << "\t\tregroup . . ." << endl;
      lout << "\t\tregroup . . ." << endl;

      time_t start_03 = time(NULL);
      double local_cost = regroup(local, lout);
      time_t finish_03 = time(NULL);

      cerr << "\t\ttime = " << finish_03 - start_03 << " seconds" << endl;
      lout << "\t\ttime = " << finish_03 - start_03 << " seconds" << endl;

      if(local_cost < best_cost)
      {
        best = local;
        best_cost = local_cost;
        cerr << "\t\taccepted" << endl;
        lout << "\t\taccepted" << endl;
      }
      else
      {
        cerr << "\t\trejected" << endl;
        lout << "\t\trejected" << endl;
      }
    }

    time_t finish_01 = time(NULL);

    global = best;
    new_cost = best_cost;

    cerr << "\tnew cost = " << new_cost << endl;
    lout << "\tnew cost = " << new_cost << endl;

    cerr << "\ttime = " << finish_01 - start_01 << " seconds" << endl;
    lout << "\ttime = " << finish_01 - start_01 << " seconds" << endl;
  }

  return global;
}

Multiclustering manual_search(Data & data, ostream & lout)
// Library facilities used: none
{
  Multiclustering local = Multiclustering(&data, &lout);
  Indexer::tuple_t dimension(data.ways()); 
  dimension[0] = -1; dimension[1] = -1;
  local.print_2D_slice(dimension, cout); cerr << local.cost() << endl;
  local.add_cluster(0);
  local.print_2D_slice(dimension, cout); cerr << local.cost() << endl;
  local.add_cluster(0);
  local.print_2D_slice(dimension, cout); cerr << local.cost() << endl;
  local.add_cluster(1);
  local.print_2D_slice(dimension, cout); cerr << local.cost() << endl;
  return local;
}

int main(int argc, char ** argv)
{
  // Check command-line arguments
  if(argc != 2) {cerr << "usage: " << argv[0] << " dir" << endl; exit(1);}

  // Set up prgoram name
  string prog_name("multi-clustering");
  cout << "Running " << prog_name << endl << endl;

  // Set up input_dir
  string input_dir(argv[1]);
  cout << "input " << input_dir << endl << endl;

  // Create output directory
  string output_dir = create_output_dir(input_dir, prog_name);
  cout << "output " << output_dir << endl << endl;

  // Open log file
  cout << "opening log file " << LOG_FILE << " . . . ";
  string log_pathname = string(output_dir + LOG_FILE);
  ofstream lout(log_pathname.c_str());
  if(lout.fail()) {cout << "error opening " << log_pathname << endl; exit(1);}
  cout << "done" << endl << endl;

  // Load data
  Data data(input_dir, DATA_FILE, LABELS_FILE);
  cout << "loading data . . . ";
  lout << "loading data . . . ";
  time_t start = time(NULL);
  data.load();
  time_t finish = time(NULL);
  cout << bytes(data.matrix.size() * sizeof(Data::T)) << " matrix, "
    << bytes(data.labels.size() * sizeof(string)) << " labels ";
  lout << bytes(data.matrix.size() * sizeof(Data::T)) << " matrix, "
    << bytes(data.labels.size() * sizeof(string)) << " labels ";
  cout << "done " << finish - start << " seconds" << endl << endl;
  lout << "done " << finish - start << " seconds" << endl << endl;

  // define plane to print
  Indexer::tuple_t dimension(data.ways(), 0);
  dimension[0] = -1;
  dimension[1] = -1;

  // Print matrix
  data.matrix.print_2D_slice(dimension, cout); cout << endl;
  data.matrix.print_2D_slice(dimension, output_dir + MATRIX_FILE);

  //vector<Data::T> symmetric_mask = symmetric_encode_data(data);
  //symmetric_encode_UN_world_trade_data(data);

  //// Print matrix
  //data.matrix.print_2D_slice(dimension, cout); cout << endl;
  //data.matrix.print_2D_slice(dimension, output_dir + MATRIX_FILE);

  //cerr.precision(numeric_limits<double>::digits10 + 3);
  cerr << "crossassociation search . . ." << endl;
  lout << "crossassociation search . . ." << endl;

  start = time(NULL);
  Multiclustering solution = crossassociation_search(data, lout);
  finish = time(NULL);

  cerr << finish - start << " seconds" << endl;

  //Multiclustering solution = manual_search(data, lout);

  //// get original data back
  //for(int i = 0; i < int(symmetric_mask.size()); ++i)
  //  if(symmetric_mask[i]) data.matrix[i] = 0;

  // print
  cerr << endl << "solution . . ." << endl;
  solution.print_2D_slice(dimension, cout);
  cerr << solution.cost() << endl;
  solution.print_blocked_matrix_2D(string(output_dir + BLOCKED_MATRIX_FILE));
  solution.print_model_2D(dimension, cout);
  solution.print_model_2D(dimension, string(output_dir + BLOCK_MODEL_FILE));
  solution.print_clusterings(output_dir);
  solution.print_block_densities(string(output_dir + BLOCK_DENSITIES_FILE));
  lout << "cost = " << solution.cost() << endl;
  lout << "time = " << finish - start << " seconds" << endl;
  lout.close();

  return 0;
}
