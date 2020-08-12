/**
 * \copyright 2011 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SPATIALHASHMAP_HH_
#define SPATIALHASHMAP_HH_

#include <functional>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "../Core/Definitions.hh"
#include "ThreadUtils.hh"

//#include "tbb/pipeline.h"
//#include "tbb/mutex.h"

namespace strandsim {

template <typename DataT, typename IteratorT, bool GenFromAABB>
class SpatialHashMap {
 private:
  inline static double DEFAULT_CELL_SIZE() { return 1.; }
  inline static double MAX_MEAN_OBJECTS_PET_CELL() { return 2.; }
  const static unsigned MAX_CELLS_PER_OBJECT = 1000;

 public:
  typedef SpatialHashMap<DataT, IteratorT, GenFromAABB> SpatialHashMapT;

  //! Sample struct : <One colliding object, one subsample ID (eg segment)>
  struct Sample : public std::pair<DataT *, IteratorT> {
   private:
    std::hash<DataT *> m_ptrHasher;
    std::hash<IteratorT> m_iterHasher;

   public:
    long hash() const {
      return m_ptrHasher(this->first) ^ m_iterHasher(this->second);
    }
  };

  // Map First object -> Second object -> Pair of vertices
  class Result;

  explicit SpatialHashMap(const Scalar cellSize = DEFAULT_CELL_SIZE())
      : m_minCellSize(cellSize),
        m_cellSize(cellSize),
        m_invCellSize(1. / cellSize) {}

  ~SpatialHashMap() {}

  void setCellSize(const double cellSize) {
    m_minCellSize = cellSize;
    m_cellSize = m_minCellSize;
    m_invCellSize = 1. / m_cellSize;
    clear();
  }

  void setMinCellSize(const double minCellSize) {
    m_minCellSize = minCellSize;
    if (m_minCellSize > 1.1 * m_cellSize) {
      m_cellSize = m_minCellSize;
      m_invCellSize = 1. / m_cellSize;
      clear();
    }
  }

  //! Adds a single object, sampling from begin to end
  void addObject(DataT &object, const IteratorT &begin, const IteratorT &end);
  //! Removes a single object
  void removeObject(DataT &object);
  //! Updates a single object, sampling from begin to end
  void updateObject(DataT &object, const IteratorT &begin,
                    const IteratorT &end);
  //! Batch updates a list of objects, sampling from DataT::subsamples_begin()
  //! to DataT::subsamples_end()
  // with | DataT::subsamples_begin() - DataT::subsamples_end() | <=
  // maxSubsamples
  void batchUpdate(const std::vector<DataT *> &object,
                   const size_t maxSubSamples);

  //! Compute possible collision pairs on a list of objects (or all if objects
  //! is NULL )
  /*! \param selfCollisions Whether to include self-collisions as well */
  void compute(Result &result, const std::list<DataT *> *objects = NULL) const;

  //! Clears the hashmap
  void clear();

  size_t size() const { return m_hashMap.size(); }

  // Footrint : Non empty -cells
  //! Adds a single object to the footprint, sampling from begin to end
  void addToFootPrint(DataT &object, const IteratorT &begin,
                      const IteratorT &end);
  //! Check if an object is at least partially in the footprint
  bool isInFootPrint(DataT &object, const IteratorT &begin,
                     const IteratorT &end) const;

  class Cell {
   private:
    int m_i;
    int m_j;
    int m_k;

    long m_hash;

   public:
    Cell(int i, int j, int k)
        : m_i(i),
          m_j(j),
          m_k(k),
          m_hash((m_i * 73856093) ^ (m_j * 19349663) ^ (m_k * 83492791)) {}

    bool operator==(const Cell &rhs) const {
      return m_i == rhs.m_i && m_j == rhs.m_j && m_k == rhs.m_k;
    }

    long hash() const { return m_hash; }

    struct Hasher {
      long operator()(const Cell &k) const { return k.hash(); }
    };

    struct Comparer {
      bool operator()(const Cell &lhs, const Cell &rhs) const {
        return lhs == rhs;
      }
    };

    enum GeneratorType { GENERATOR = GenFromAABB };
  };

  typedef std::unordered_set<Cell, typename Cell::Hasher,
                             typename Cell::Comparer>
      HashSet;
  typedef std::map<DataT *, std::vector<IteratorT> > Samples;
  typedef std::unordered_map<Cell, Samples, typename Cell::Hasher,
                             typename Cell::Comparer>
      HashMap;
  typedef std::unordered_map<DataT *, HashSet> ObjectCells;

 private:
  Scalar m_minCellSize;
  Scalar m_cellSize;
  Scalar m_invCellSize;

  ObjectCells m_objectCells;

  HashMap m_hashMap;
  HashSet m_footPrint;

 public:
  // Parallelisation of rasterization using tbb::pipeline formalism
  class ParallelAdder_Buffer;     //!< Buffer for chunks of objects
  class ParallelAdder_Input;      //!< Reads input, allocates buffers
  class ParallelAdder_Transform;  //!< Proper rasterizer
  class ParallelAdder_Output;
  //!< Writes cells to hashmap, delete buffers

  class Result {
   public:
    typedef std::unordered_map<
        DataT *, std::unordered_map<
                     DataT *, std::set<std::pair<IteratorT, IteratorT> > > >
        Collisions;

    // Public API
    Collisions &collisions() {
      if (!done()) next(m_collisions);
      return m_collisions;
    }

    bool done() const { return m_cur == m_objectsIters.size(); }
    bool next(Collisions &collisions);

    Result(bool doSelfCollisions = false, unsigned batchSize = 0)
        : m_cur(0),
          m_batchSize(batchSize),
          m_doSelfCollisions(doSelfCollisions),
          m_checkAllObjectsInCells(false),
          m_hashMap(NULL) {}

   private:
    void reset(const SpatialHashMap *hashMap,
               std::vector<typename ObjectCells::const_iterator> &objectsIters,
               bool checkAllObjectsInCells) {
      m_hashMap = hashMap;
      m_cur = 0;
      m_objectsIters.swap(objectsIters);
      m_checkAllObjectsInCells = checkAllObjectsInCells;
    }

    unsigned m_cur;
    unsigned m_batchSize;
    bool m_doSelfCollisions;
    bool m_checkAllObjectsInCells;  // If true, dont' use symmetry property
    const SpatialHashMap *m_hashMap;

    Collisions m_collisions;
    std::vector<typename ObjectCells::const_iterator> m_objectsIters;

    friend class SpatialHashMap;
  };

  void compute(const typename ObjectCells::const_iterator objectCellsIt,
               bool selfCollisions, bool checkAllObjectsInCells,
               typename Result::Collisions &collisions) const;
  void compute(typename HashMap::const_iterator,
               typename Samples::const_iterator object1, bool selfCollisions,
               bool checkAllObjectsInCells,
               typename Result::Collisions &collisions) const;

  //! Generates cells from a 3D segment
  template <typename ContainerT>
  static void generateCellsFromSegment(const DataT &object, const IteratorT &it,
                                       const Scalar invCellSize,
                                       ContainerT &cells);

  //! Generates cells from a 3d axis-aligned bounding box
  template <typename ContainerT>
  static void generateCellsFromAABB(DataT &object, const IteratorT &it,
                                    const Scalar invCellSize,
                                    ContainerT &cells);

  friend class ParallelAdder_Buffer;     //!< Buffer for chunks of objects
  friend class ParallelAdder_Input;      //!< Reads input, allocates buffers
  friend class ParallelAdder_Transform;  //!< Proper rasterizer
  friend class ParallelAdder_Output;
  friend class Result;
};

template <typename DataT, typename IteratorT, typename ContainerT>
void generateHashMapCells(DataT &object, const IteratorT &it,
                          const Scalar invCellSize, ContainerT &cells,
                          const typename SpatialHashMap<
                              DataT, IteratorT, false>::Cell::GeneratorType) {
  SpatialHashMap<DataT, IteratorT, false>::generateCellsFromSegment(
      object, it, invCellSize, cells);
}

template <typename DataT, typename IteratorT, typename ContainerT>
void generateHashMapCells(DataT &object, const IteratorT &it,
                          const Scalar invCellSize, ContainerT &cells,
                          const typename SpatialHashMap<
                              DataT, IteratorT, true>::Cell::GeneratorType) {
  SpatialHashMap<DataT, IteratorT, true>::generateCellsFromAABB(
      object, it, invCellSize, cells);
}

template <typename DataT, typename IteratorT, typename ContainerT>
void generateHashMapCells(DataT &object, const IteratorT &it,
                          const Scalar invCellSize, ContainerT &cells) {
  typedef typename ContainerT::value_type CellT;
  generateHashMapCells<DataT, IteratorT, ContainerT>(object, it, invCellSize,
                                                     cells, CellT::GENERATOR);
}

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::clear() {
  m_objectCells.clear();
  m_hashMap.clear();
  m_footPrint.clear();
}

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::addObject(DataT &object,
                                                        const IteratorT &begin,
                                                        const IteratorT &end) {
  HashSet objectCells;

  std::vector<std::vector<Cell> > sampleCells(end - begin);

  for (IteratorT it = begin; it < end; ++it) {
    generateHashMapCells<DataT, IteratorT>(object, it, m_invCellSize,
                                           sampleCells[it]);
  }

  for (IteratorT it = begin; it < end; ++it) {
    // Concatenate sample's cells to object's cells
    objectCells.insert(sampleCells[it].begin(), sampleCells[it].end());

    // Add sample in each cell of spatial hash map
    for (typename std::vector<Cell>::const_iterator cell =
             sampleCells[it].begin();
         cell != sampleCells[it].end(); ++cell) {
      m_hashMap[*cell][&object].push_back(it);
    }
  }

  m_objectCells[&object].swap(objectCells);

  // std::cout << "Add -> Current hash map size:" << m_hashMap.size() <<
  // std::endl ;
}

// Footrint : Non empty -cells
//! Adds a single object to the footprint, sampling from begin to end
template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::addToFootPrint(
    DataT &object, const IteratorT &begin, const IteratorT &end) {
  std::vector<std::list<Cell> > sampleCells(end - begin);

  for (IteratorT it = begin; it < end; ++it) {
    generateHashMapCells<DataT, IteratorT>(object, it, m_invCellSize,
                                           sampleCells[it - begin]);
  }

  for (IteratorT it = begin; it < end; ++it) {
    // Concatenate sample's cells to object's cells
    m_footPrint.insert(sampleCells[it - begin].begin(),
                       sampleCells[it - begin].end());
  }
}

template <typename DataT, typename IteratorT, bool GAABB>
bool SpatialHashMap<DataT, IteratorT, GAABB>::isInFootPrint(
    DataT &object, const IteratorT &begin, const IteratorT &end) const {
  std::vector<Cell> sampleCells;

  for (IteratorT it = begin; it < end; ++it) {
    sampleCells.clear();
    generateHashMapCells<DataT, IteratorT>(object, it, m_invCellSize,
                                           sampleCells);

    for (typename std::vector<Cell>::const_iterator cell = sampleCells.begin();
         cell != sampleCells.end(); ++cell) {
      if (m_footPrint.find(*cell) != m_footPrint.end()) {
        return true;
      }
    }
  }

  return false;
}

template <typename DataT, typename IteratorT, bool GAABB>
class SpatialHashMap<DataT, IteratorT, GAABB>::ParallelAdder_Buffer {
 public:
  unsigned m_begin;
  unsigned m_end;

  bool m_free;

  // Set of object's cells
  std::vector<HashSet> m_objectCells;

  typedef std::unordered_map<Cell, std::vector<IteratorT>,
                             typename Cell::Hasher, typename Cell::Comparer>
      HashMap;

  // Object's hash map : list of all subsamples inside each cell
  // actually m_objectCells should be the keys of m_objectMaps, bu can't find
  // a way to specify that using unordered_map
  std::vector<HashMap> m_objectMaps;

  ParallelAdder_Buffer() {}

  ~ParallelAdder_Buffer() {}

  void assign(unsigned begin, unsigned end) {
    m_begin = begin;
    m_end = end;
    m_free = false;
  }
  void release() {
    clear();
    m_free = true;
  }
  bool isFree() const { return m_free; }

 private:
  void clear() {
    for (int i = 0; i < m_objectMaps.size(); ++i) {
      m_objectMaps[i].clear();
      m_objectCells[i].clear();
    }
  }
};

template <typename DataT, typename IteratorT, bool GAABB>
class SpatialHashMap<DataT, IteratorT,
                     GAABB>::ParallelAdder_Input  //: public tbb::filter
{
 private:
  const std::vector<DataT *> &m_objects;

  unsigned m_curIndex;
  unsigned m_chunkSize;

  std::list<ParallelAdder_Buffer *> m_buffers;

 public:
  ParallelAdder_Input(const std::vector<DataT *> &objects, unsigned chunkSize)
      :  // tbb::filter(true),
        m_objects(objects),
        m_curIndex(0),
        m_chunkSize(chunkSize) {}
  virtual ~ParallelAdder_Input() {
    //        std::cout << "Did everythin with " << m_buffers.size() << "
    //        buffers. Great. " << std::endl ;
    for (auto bufIt = m_buffers.begin(); bufIt != m_buffers.end(); ++bufIt) {
      delete *bufIt;
    }
  }

  void *operator()(void *data) {
    ParallelAdder_Buffer *b = NULL;
    if (m_curIndex < m_objects.size()) {
      unsigned next =
          std::min(m_curIndex + m_chunkSize, (unsigned)m_objects.size());

      for (auto bufIt = m_buffers.begin(); bufIt != m_buffers.end(); ++bufIt) {
        if ((*bufIt)->isFree()) {
          b = *bufIt;
          break;
        }
      }
      if (!b) {
        b = new ParallelAdder_Buffer();
        m_buffers.push_back(b);
      }

      b->assign(m_curIndex, next);
      m_curIndex = next;
    }

    return b;
  }
};

template <typename DataT, typename IteratorT, bool GAABB>
class SpatialHashMap<DataT, IteratorT,
                     GAABB>::ParallelAdder_Transform  //: public tbb::filter
{
 private:
  const std::vector<DataT *> &m_objects;

  const Scalar m_invCellSize;
  const Scalar m_maxSubSamples;

 public:
  ParallelAdder_Transform(const std::vector<DataT *> &objects,
                          const size_t maxSubsamples,
                          const Scalar invCellSize)
      :  // tbb::filter(false),
        m_objects(objects),
        m_invCellSize(invCellSize),
        m_maxSubSamples(maxSubsamples) {}
  virtual ~ParallelAdder_Transform() {}

  void *operator()(void *data) {
    ParallelAdder_Buffer &b = *static_cast<ParallelAdder_Buffer *>(data);

    // std::cout << "Start " << b.m_begin << " : " << b.m_end << "\n" ;

    b.m_objectCells.resize(b.m_end - b.m_begin);
    b.m_objectMaps.resize(b.m_end - b.m_begin);

    std::vector<std::vector<Cell> > sampleCells(m_maxSubSamples);

    for (int i = b.m_begin; i < b.m_end; ++i) {
      const IteratorT samplesBegin = m_objects[i]->subsamples_begin();
      const IteratorT samplesEnd = m_objects[i]->subsamples_end();

      const unsigned curIdx = i - b.m_begin;

      for (IteratorT it = samplesBegin; it < samplesEnd; ++it) {
        sampleCells[it].clear();
        generateHashMapCells<DataT, IteratorT>(*m_objects[i], it, m_invCellSize,
                                               sampleCells[it]);

        // Concatenate sample's cells to object's cells
        b.m_objectCells[curIdx].insert(sampleCells[it].begin(),
                                       sampleCells[it].end());
      }

      // Fill object's map
      for (IteratorT it = samplesBegin; it < samplesEnd; ++it) {
        for (typename std::vector<Cell>::const_iterator cell =
                 sampleCells[it].begin();
             cell != sampleCells[it].end(); ++cell) {
          b.m_objectMaps[curIdx][*cell].push_back(it);
        }
      }
    }
    // std::cout << "Done " << b.m_begin << " : " << b.m_end << "\n" ;

    return &b;
  }

  class StaticArg {
   public:
    StaticArg() : m_object(NULL), m_data(NULL) {}

    void setFunctor(ParallelAdder_Transform &f) { m_object = &f; }

    void setData(void *data) { m_data = data; }
    void *data() { return m_data; }

    void operator()() { (*this->m_object)(this->m_data); }

    static void *call(void *packedData) {
      StaticArg *staticArg = static_cast<StaticArg *>(packedData);
      return (*staticArg->m_object)(staticArg->m_data);
    }

   private:
    ParallelAdder_Transform *m_object;
    void *m_data;
  };
};

template <typename DataT, typename IteratorT, bool GAABB>
class SpatialHashMap<DataT, IteratorT,
                     GAABB>::ParallelAdder_Output  //: public tbb::filter
{
 private:
  SpatialHashMap<DataT, IteratorT, GAABB> &m_hashMap;
  const std::vector<DataT *> &m_objects;

  bool m_doCleanUp;

  // tbb::mutex m_mutex ;

 public:
  ParallelAdder_Output(SpatialHashMap<DataT, IteratorT, GAABB> &hashMap,
                       const std::vector<DataT *> &objects, bool doCleanUp)
      :  // tbb::filter( true ),
        m_hashMap(hashMap),
        m_objects(objects),
        m_doCleanUp(doCleanUp) {}
  virtual ~ParallelAdder_Output() {}

  void *operator()(void *data) {
    ParallelAdder_Buffer &b = *static_cast<ParallelAdder_Buffer *>(data);

    // std::cout << "Start writing " << b.m_begin << " : " << b.m_end << "\n" ;
#pragma omp critical
    {
      // typename tbb::mutex::scoped_lock lock(m_mutex) ;

      for (int i = b.m_begin; i < b.m_end; ++i) {
        const unsigned curIdx = i - b.m_begin;

        // Fill global hashmap from individual objects' maps
        for (typename ParallelAdder_Buffer::HashMap::iterator iter =
                 b.m_objectMaps[curIdx].begin();
             iter != b.m_objectMaps[curIdx].end(); ++iter) {
          m_hashMap.m_hashMap[iter->first][m_objects[i]].swap(iter->second);
        }

        // If the hashmap was not empty, we have to remove objects from cells in
        // which they are not in anymore
        HashSet &objectCells = m_hashMap.m_objectCells[m_objects[i]];
        if (m_doCleanUp) {
          for (typename HashSet::const_iterator cell = objectCells.begin();
               cell != objectCells.end(); ++cell) {
            if (b.m_objectCells[curIdx].end() ==
                b.m_objectCells[curIdx].find(*cell)) {
              m_hashMap.m_hashMap.find(*cell)->second.erase(m_objects[i]);
            }
          }
        }
        objectCells.swap(b.m_objectCells[curIdx]);
      }
    }

    b.release();
    return NULL;
  }
};

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::batchUpdate(
    const std::vector<DataT *> &objects, const size_t maxSubSamples) {
  const size_t nObj = objects.size();
  const bool doCleanUp = m_hashMap.size() > 0;

  {
    ParallelAdder_Input input(objects, 10);
    ParallelAdder_Output output(*this, objects, doCleanUp);
    ParallelAdder_Transform transform(objects, maxSubSamples, m_invCellSize);

    void *data = NULL;
    // Sequential / OMP version
#pragma omp parallel private(data)
    do {
#pragma omp critical
      data = input(NULL);
    } while (data && !output(transform(data)));

    // tbb version -- not worth it
    //    tbb::pipeline pipeline ;

    //     pipeline.add_filter( input ) ;

    //     pipeline.add_filter( transform ) ;

    //     pipeline.add_filter( output ) ;

    //     pipeline.run( 8 ) ;
    //     pipeline.clear() ;
  }

  // Clean up: remove empty cells
  if (doCleanUp) {
    for (typename HashMap::iterator iter = m_hashMap.begin();
         iter != m_hashMap.end(); ++iter) {
      if (iter->second.size() == 0) {
        m_hashMap.erase(iter);
      }
    }
  }

  //    std::cout << "batchUpdate " << nObj << "-> Current hash map size:" <<
  //    m_hashMap.size() << ", cell size: " << m_cellSize << std::endl ;
  //    std::cout << m_hashMap.load_factor() << std::endl ;

  // Check if over-crowded
  // -> Decrease cell size if it is ( should staty higher than cell thickness )

  if (m_cellSize >= 2 * m_minCellSize) {
    unsigned objSum = 0;
    for (typename HashMap::const_iterator iter = m_hashMap.begin();
         iter != m_hashMap.end(); ++iter) {
      objSum += iter->second.size();
    }
    const double meanOPC = ((double)objSum) / m_hashMap.size();

    // std::cout << "Cell size: " << m_cellSize << "; Mean objects per cell: "
    // << meanOPC << std::endl ;

    if (meanOPC > MAX_MEAN_OBJECTS_PET_CELL()) {
      m_cellSize = std::max(m_minCellSize, m_cellSize * (Scalar).5);
      m_invCellSize = 1. / m_cellSize;

      clear();
      batchUpdate(objects, maxSubSamples);
    }
  }
}

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::removeObject(DataT &object) {
  HashSet &objectCells = m_objectCells[&object];

  // std::cout << "Remove -> Num of Object's cells:" << objectCells.size() <<
  // std::endl ;

  for (typename HashSet::iterator cell = objectCells.begin();
       cell != objectCells.end(); ++cell) {
    typename HashMap::iterator mapCell = m_hashMap.find(*cell);
    if (mapCell != m_hashMap.end()) {
      Samples &samples = mapCell->second;

      if (samples.size() == 1) {
        m_hashMap.erase(*cell);
      } else {
        samples.erase(&object);
      }
    }
  }
  // std::cout << "Remove -> Current hash map size:" << m_hashMap.size() <<
  // std::endl ;

  objectCells.clear();
}

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::updateObject(
    DataT &object, const IteratorT &begin, const IteratorT &end) {
  removeObject(object);
  addObject(object, begin, end);
}

template <typename DataT, typename IteratorT, bool GAABB>
bool SpatialHashMap<DataT, IteratorT, GAABB>::Result::next(
    Collisions &collisions) {
  if (!m_hashMap) return false;

  unsigned begin, end;

#pragma omp critical(nextCollisionBatch)
  {
    begin = m_cur;
    end = m_batchSize
              ? std::min(m_cur + m_batchSize, (unsigned)m_objectsIters.size())
              : (unsigned)m_objectsIters.size();
    m_cur = end;
  }

  if (begin == end) return false;

  collisions.clear();
  for (int iterIdx = begin; iterIdx < end; ++iterIdx) {
    collisions[m_objectsIters[iterIdx]->first];  // Create external hashmap
                                                 // structure
  }

#pragma omp parallel for
  for (int iterIdx = begin; iterIdx < end; ++iterIdx) {
    m_hashMap->compute(m_objectsIters[iterIdx], m_doSelfCollisions,
                       m_checkAllObjectsInCells, collisions);
  }

  return true;
}

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::compute(
    const typename ObjectCells::const_iterator objectCellsIt,
    bool selfCollisions, bool checkAllObjectsInCells,
    typename Result::Collisions &collisions) const {
  const HashSet &cells = objectCellsIt->second;
  for (typename HashSet::const_iterator cell = cells.begin();
       cell != cells.end(); ++cell) {
    typename HashMap::const_iterator mapCell = m_hashMap.find(*cell);
    if (mapCell != m_hashMap.end()) {
      const Samples &samples = mapCell->second;
      const typename Samples::const_iterator object1 =
          samples.find(objectCellsIt->first);
      if (object1 != samples.end()) {
        compute(mapCell, object1, selfCollisions, checkAllObjectsInCells,
                collisions);
      }
    }
  }
}

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::compute(
    typename HashMap::const_iterator iter,
    typename Samples::const_iterator object1, bool selfCollisions,
    bool checkAllObjectsInCells,
    typename Result::Collisions &collisions) const {
  const Samples &samples = iter->second;

  auto &first = collisions[object1->first];

  const typename Samples::const_iterator begin =
      checkAllObjectsInCells ? samples.begin() : object1;

  for (typename Samples::const_iterator object2 = begin;
       object2 != samples.end(); ++object2) {
    if (!selfCollisions && object2 == object1) continue;

    auto &second = first[object2->first];

    for (int index1 = 0; index1 < object1->second.size(); ++index1) {
      const IteratorT elem1 = object1->second[index1];

      unsigned startIndex = (object2 == object1) ? index1 + 2 : 0;

      for (int index2 = startIndex; index2 < object2->second.size(); ++index2) {
        second.insert(std::make_pair(elem1, object2->second[index2]));
      }
    }
  }
}

template <typename DataT, typename IteratorT, bool GAABB>
void SpatialHashMap<DataT, IteratorT, GAABB>::compute(
    Result &result, const std::list<DataT *> *objects) const {
  std::vector<typename ObjectCells::const_iterator> objectsIters;

  if (objects) {
    objectsIters.resize(objects->size());

    int i = 0;
    for (typename std::list<DataT *>::const_iterator object = objects->begin();
         object != objects->end(); ++object) {
      typename ObjectCells::const_iterator objectCellsIt =
          m_objectCells.find(*object);

      objectsIters[i++] = objectCellsIt;
    }

  } else {
    const unsigned nObjects(m_objectCells.size());
    objectsIters.resize(nObjects);

    int i = 0;
    for (auto objectCellsIt = m_objectCells.begin();
         objectCellsIt != m_objectCells.end(); ++objectCellsIt) {
      objectsIters[i++] = objectCellsIt;
    }
  }

  result.reset(this, objectsIters, objects);
}

template <typename DataT, typename IteratorT, bool GAABB>
template <typename ContainerT>
void SpatialHashMap<DataT, IteratorT, GAABB>::generateCellsFromSegment(
    const DataT &object, const IteratorT &it, const Scalar invCellSize,
    ContainerT &cells) {
  Vec3x start, end;
  object.getSegment(it, start, end);

  const Vec3x dP = (end - start) * invCellSize;

  unsigned kMax = 0;
  Scalar dMax = std::fabs(dP[0]);

  for (unsigned k = 1; k < 3; ++k) {
    const Scalar a = std::fabs(dP[k]);
    if (a > dMax) {
      kMax = k;
      dMax = a;
    }
  }

  const int sign = dP[kMax] > 0 ? 1 : -1;
  const unsigned n = std::ceil(dMax) + 1;

  if (n > MAX_CELLS_PER_OBJECT) {
    std::cerr << "SHM - Warning: Unusually big object " << std::endl;
    return;
  }

  const Vec3x dir = dP / dMax;

  Vec3x cur = start * invCellSize;

  uint16_t ip = kMax;
  uint16_t jp = (kMax + 1) % 3;
  uint16_t kp = (kMax + 2) % 3;

  int permuted[3] = {
      (int)std::floor(cur[ip]),
      (int)std::floor(cur[jp]),
      (int)std::floor(cur[kp]),
  };

  int &rip = permuted[((3 - ip) % 3)];
  int &rjp = permuted[((3 - jp) % 3)];
  int &rkp = permuted[((3 - kp) % 3)];

  for (int i = 0; i < n; ++i) {
    cells.push_back(Cell(rip, rjp, rkp));

    cur += dir;

    int jip1 = std::floor(cur[jp]);
    int kip1 = std::floor(cur[kp]);

    // One of the other directions have changed cell ; add cell at altitude i
    if (jip1 != permuted[1] || kip1 != permuted[2]) {
      // Both other directions have changed cells ; add yet another cell
      // For now add two, but FIXME !
      if (jip1 != permuted[1] && kip1 != permuted[2]) {
        const int oldj = jip1;
        permuted[1] = jip1;
        cells.push_back(Cell(rip, rjp, rkp));

        permuted[1] = oldj;
        permuted[2] = kip1;
        cells.push_back(Cell(rip, rjp, rkp));
      }

      permuted[1] = jip1;
      permuted[2] = kip1;

      cells.push_back(Cell(rip, rjp, rkp));
    }

    permuted[0] += sign;
  }
}

template <typename DataT, typename IteratorT, bool GAABB>
template <typename ContainerT>
void SpatialHashMap<DataT, IteratorT, GAABB>::generateCellsFromAABB(
    DataT &object, const IteratorT &it, const Scalar invCellSize,
    ContainerT &cells) {
  Vec3x min, max;
  object.getAABB(it, min, max);

  int imin = ((int)std::ceil(min[0] * invCellSize)) - 1;
  int jmin = ((int)std::ceil(min[1] * invCellSize)) - 1;
  int kmin = ((int)std::ceil(min[2] * invCellSize)) - 1;
  int imax = ((int)std::floor(max[0] * invCellSize)) + 1;
  int jmax = ((int)std::floor(max[1] * invCellSize)) + 1;
  int kmax = ((int)std::floor(max[2] * invCellSize)) + 1;

  // Safeguard for diverging objects
  //    if ( ( imax - imin ) * ( jmax - jmin ) * ( kmax - kmin ) >
  //    MAX_CELLS_PER_OBJECT )
  //    {
  //        std::cerr << "SHM - Warning: Unusually big object " << std::endl;
  //        return;
  //    }

  //    cells.reserve( ( imax - imin ) * ( jmax - jmin ) * ( kmax - kmin ) ) ;

  for (int i = imin; i < imax; ++i) {
    for (int j = jmin; j < jmax; ++j) {
      for (int k = kmin; k < kmax; ++k) {
        cells.push_back(Cell(i, j, k));
      }
    }
  }
}

}  // namespace strandsim

#endif /* SPATIALHASHMAP_HH_ */
