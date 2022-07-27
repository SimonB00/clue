#ifndef CLUEAlgo_h
#define CLUEAlgo_h

// C/C++ headers
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <chrono>

#include "LayerTiles.h"
#include "Points.h"

template <uint8_t N>
class CLUEAlgo{
public:
  CLUEAlgo(float dc, float rhoc, float outlierDeltaFactor, bool verbose) {
    dc_ = dc; 
    rhoc_ = rhoc;
    outlierDeltaFactor_ = outlierDeltaFactor;
    verbose_ = verbose;
  }
  ~CLUEAlgo(){} 
    
  // public variables
  float dc_;  // cut-off distance in the calculation of local density
  float rhoc_;  // minimum density to promote a point as a seed or the maximum density to demote a point as an outlier
  float outlierDeltaFactor_;
  bool verbose_;

  float dm = outlierDeltaFactor_ * dc_; // separation requirement for seeds (I suppose)
    
  Points<N> points_;
  
  //bool setPoints(int n, float* x, float* y, int* layer, float* weight) {
  bool setPoints(int n, std::array<std::vector<float>,N>& coordinates, int* layer, float* weight) {
    
    points_.clear();
    // input variables
    for(int i = 0; i < n; ++i) {
	    //points_.x.push_back(x[i]);
	    //points_.y.push_back(y[i]);
	    for(int j = 0; j != N; ++j) {
        points_.coordinates_[j].push_back(coordinates[j][i]);
      }
      points_.layer.push_back(layer[i]);
	    points_.weight.push_back(weight[i]);
    }

    points_.n = points_.coordinates_[0].size();
    if(points_.n == 0)
      return 1;

    // result variables
    points_.rho.resize(points_.n,0);
    points_.delta.resize(points_.n,std::numeric_limits<float>::max());
    points_.nearestHigher.resize(points_.n,-1);
    points_.followers.resize(points_.n);
    points_.clusterIndex.resize(points_.n,-1);
    points_.isSeed.resize(points_.n,0);
    return 0;
  }

  void clearPoints(){ points_.clear(); }

  void makeClusters() {
    std::array<LayerTiles<N>, NLAYERS> allLayerTiles;
    // start clustering
    auto start = std::chrono::high_resolution_clock::now();

    prepareDataStructures(allLayerTiles);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "--- prepareDataStructures:     " << elapsed.count() *1000 << " ms\n";

    start = std::chrono::high_resolution_clock::now();
    calculateLocalDensity(allLayerTiles);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "--- calculateLocalDensity:     " << elapsed.count() *1000 << " ms\n";

    start = std::chrono::high_resolution_clock::now();
    calculateDistanceToHigher(allLayerTiles);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "--- calculateDistanceToHigher: " << elapsed.count() *1000 << " ms\n";

    findAndAssignClusters();
  };

  void infoSeeds();
  void infoHits();

  void for_recursion(int N_, std::vector<int> &base_vector,  std::vector<int> &dim_min, std::vector<int> &dim_max, LayerTiles<N>& lt_, int point_id) {
    if(!N) {
        int binId = lt_.getGlobalBinByBin(base_vector);
        // get the size of this bin
        int binSize = lt_[binId].size();
        
        // iterate inside this bin
        for (int binIter = 0; binIter < binSize; ++binIter) {
          int j = lt_[binId][binIter];
          // query N_{dc_}(i)
          float dist_ij = distance(point_id, j);
          if(dist_ij <= dc_) {
            // sum weights within N_{dc_}(i)
            points_.rho[point_id] += (point_id == j ? 1.f : 0.5f) * points_.weight[j];
          }
        } // end of interate inside this bin
    }
    for(int i = dim_min[dim_min.size() - N]; i < dim_max[dim_max.size() - N]; ++i) {
        base_vector[base_vector.size()-N] = i;
        for_recursion(N-1, base_vector, dim_min, dim_max, lt_, point_id);
    }
  };

  //for_recursion used for the function calculateDistanceToHigher
  void for_recursion_DistanceToHigher(int N_, std::vector<int> &base_vector,  std::vector<int> &dim_min, std::vector<int> &dim_max, 
    LayerTiles<N>& lt_, float rho_i, float delta_i, int nearestHigher_i, int point_id) {
      if(!N){
          int binId = lt_.getGlobalBinByBin(base_vector);
          // get the size of this bin
          int binSize = lt_[binId].size();

          // iterate inside this bin
          for (int binIter = 0; binIter < binSize; binIter++) {
            int j = lt_[binId][binIter];
            // query N'_{dm}(i)
            bool foundHigher = (points_.rho[j] > rho_i);
            // in the rare case where rho is the same, use detid
            foundHigher = foundHigher || ((points_.rho[j] == rho_i) && (j > point_id) );
            float dist_ij = distance(point_id, j);
            if(foundHigher && dist_ij <= dm) { // definition of N'_{dm}(i)
              // find the nearest point within N'_{dm}(i)
              if (dist_ij < delta_i) {
                // update delta_i and nearestHigher_i
                delta_i = dist_ij;
                nearestHigher_i = j;
              }
            }
          } // end of interate inside this bin
      }
      for(int i = dim_min[dim_min.size() - N]; i < dim_max[dim_max.size() - N]; ++i){
          base_vector[base_vector.size()-N] = i;
          for_recursion_DistanceToHigher(N-1, base_vector, dim_min, dim_max, lt_, rho_i, delta_i, nearestHigher_i, point_id);
      }
  }

  std::string getVerboseString_(unsigned it,
				float x, float y, int layer, float weight,
				float rho, float delta,
				int nh, int isseed, float clusterid,
				unsigned nVerbose) const {
    std::stringstream s;
    std::string sep = ",";
    s << it << sep << x << sep << y << sep;
    s << layer << sep << weight << sep << rho;
    if (delta <= 999)
      s << sep << delta;
    else
      s << ",999"; //convert +inf to 999 in verbose
    s << sep << nh << sep << isseed << sep << clusterid << '\n';
    return s.str();
  }
  
  void verboseResults(std::string outputFileName="cout", unsigned nVerbose=-1) const {
  //  if(verbose_)
  //    {
	//if (nVerbose==-1) nVerbose=points_.n;
  //  
	//std::string s;
	//s = "index,x,y,layer,weight,rho,delta,nh,isSeed,clusterId\n";
	//for(unsigned i=0; i<nVerbose; i++) {
	//  s += getVerboseString_(i, points_.x[i], points_.y[i], points_.layer[i],
	//			 points_.weight[i], points_.rho[i], points_.delta[i],
	//			 points_.nearestHigher[i], points_.isSeed[i],
	//			 points_.clusterIndex[i], nVerbose);
	//}
  //
	//if(outputFileName.compare("cout")==0) //verbose to screen
	//  std::cout << s << '\n;
	//else { //verbose to file
	//  std::ofstream oFile(outputFileName);
	//  oFile << s;
	//  oFile.close();
	//}
  //    }
  ;
  }
        
private:
  // private member methods
  void prepareDataStructures(std::array<LayerTiles<N>, NLAYERS> & allLayerTiles) {
    for (int i=0; i<points_.n; ++i){
      // push index of points into tiles
      std::vector<float> coords;
      for(int j = 0; j != N; ++j) {
        coords.push_back(points_.coordinates_[j][i]);
      }
      allLayerTiles[points_.layer[i]].fill(coords, i);
      // so it simply takes the layer in which the hits where detected (there is only 1 layer actually, so it should be easier),
      // divides them in tiles (bins) and saves the index of the point (hit) recorded in each of them.
    }
  };

  void calculateLocalDensity(std::array<LayerTiles<N>, NLAYERS> & allLayerTiles) {
    // loop over all points
    for(unsigned i = 0; i < points_.n; ++i) {
      LayerTiles<N>& lt = allLayerTiles[points_.layer[i]]; // there is only one layer, so this will always be the same

      // get search box
      std::array<std::vector<float>,N> minMax;
      for(int j = 0; j != N; ++j) {
        std::vector<float> partial_minMax{points_.coordinates_[j][i]-dc_,points_.coordinates_[j][i]+dc_};
        minMax[j] = partial_minMax;
      }
      std::array<int,2*N> search_box = lt.searchBox(minMax);

      // loop over bins in the search box
      std::vector<int> binVec(N);
      std::vector<int> dimMin;
      std::vector<int> dimMax;
      for(int j = 0; j != search_box.size(); ++j) {
        if(j%2 == 0) {
          dimMin.push_back(search_box[j]);
        } else {
          dimMax.push_back(search_box[j]);
        }
      }
      for_recursion(N,binVec,dimMin,dimMax,lt,i);
    } // end of loop over points
  };

  void calculateDistanceToHigher(std::array<LayerTiles<N>, NLAYERS> & allLayerTiles ) {
    // loop over all points
    for(unsigned i = 0; i < points_.n; ++i) {
      // default values of delta and nearest higher for i
      float delta_i = std::numeric_limits<float>::max();
      int nearestHigher_i = -1; // if this doesn't change, the point is either a seed or an outlier
      float rho_i = points_.rho[i];

      LayerTiles<N>& lt = allLayerTiles[points_.layer[i]];

      //get search box
      std::array<std::vector<float>,N> minMax;
      for(int j = 0; j != N; ++j) {
        std::vector<float> partial_minMax{points_.coordinates_[j][i]-dc_,points_.coordinates_[j][i]+dc_};
        minMax[j] = partial_minMax;
      }
      std::array<int,2*N> search_box = lt.searchBox(minMax);

      // loop over all bins in the search box

      std::vector<int> binVec(N);
      std::vector<int> dimMin;
      std::vector<int> dimMax;
      for(int j = 0; j != search_box.size(); ++j) {
        if(j%2 == 0) {
          dimMin.push_back(search_box[j]);
        } else {
          dimMax.push_back(search_box[j]);
        }
      }
      for_recursion_DistanceToHigher(N,binVec,dimMin,dimMax,lt, rho_i, delta_i, nearestHigher_i, i);

      points_.delta[i] = delta_i;
      points_.nearestHigher[i] = nearestHigher_i;
    } // end of loop over points
  }
  void findAndAssignClusters() {
     auto start = std::chrono::high_resolution_clock::now();

    int nClusters = 0;

    // find cluster seeds and outlier
    std::vector<int> localStack;  // this vector will contain the indexes of all the seeds
    // loop over all points
    for(unsigned i = 0; i < points_.n; i++) {
      // initialize clusterIndex
      points_.clusterIndex[i] = -1;

      float deltai = points_.delta[i];
      float rhoi = points_.rho[i];

      // determine seed or outlier 
      bool isSeed = (deltai > dc_) && (rhoi >= rhoc_);
      bool isOutlier = (deltai > outlierDeltaFactor_ * dc_) && (rhoi < rhoc_);
      if (isSeed) {
	      // set isSeed as 1
	      points_.isSeed[i] = 1;
	      // set cluster id
	      points_.clusterIndex[i] = nClusters;
	      // increment number of clusters
	      ++nClusters;
	      // add seed into local stack
	      localStack.push_back(i);
      } else if (!isOutlier) {
	      // register as follower at its nearest higher
	      points_.followers[points_.nearestHigher[i]].push_back(i);
      }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "--- findSeedAndFollowers:      " << elapsed.count() *1000 << " ms\n";

    start = std::chrono::high_resolution_clock::now();
    // expend clusters from seeds
    while (!localStack.empty()) {
      int i = localStack.back();
      auto& followers = points_.followers[i];
      localStack.pop_back();

      // loop over followers
      for(int j : followers){
        // pass id from i to a i's follower
        points_.clusterIndex[j] = points_.clusterIndex[i];
        // push this follower to localStack
        localStack.push_back(j);
      }
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "--- assignClusters:            " << elapsed.count() *1000 << " ms\n";
  };
  inline float distance(int i, int j) const {
    // 2-d distance on the layer
    if(points_.layer[i] == points_.layer[j] ) {
      //const float dx = points_.x[i] - points_.x[j];
      //const float dy = points_.y[i] - points_.y[j];
      //return std::sqrt(dx * dx + dy * dy);
      float qSum = 0.f;   // quadratic sum
      for(int k = 0; k != N; ++k) {
        qSum += std::pow(points_.coordinates_[k][i] - points_.coordinates_[k][j],2);
      }
      return std::sqrt(qSum);
    } else {
      return std::numeric_limits<float>::max();
    }
  };
};

#endif
