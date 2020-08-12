/**
 * \copyright 2012 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_MEMORY_HH
#define STRANDSIM_MEMORY_HH

#include <cstdlib>
#include <cstring>
#include <map>
#include <string>

#include "TextLog.hh"
#include "ThreadUtils.hh"

namespace strandsim {

//! Returns the current RSS from /proc/self/stat
size_t getResidentMemory();

//! Computes and prints scope-base variations of getResidentMemory().
template <typename StreamT>
class MemoryDiff {
 public:
  MemoryDiff(const std::string& name) : m_name(name), m_origRss((size_t)-1) {
    if (StreamT::would_print(g_log, "Memory")) {
      m_origRss = getResidentMemory();
    }
  }

  ~MemoryDiff() {
    if (m_origRss != ((size_t)-1)) {
      size_t newRss = getResidentMemory();
      long diff = ((long)newRss) - m_origRss;
      StreamT(g_log, "Memory") << m_name << " :: "
                               << " New: " << newRss << " Old: " << m_origRss
                               << " Diff: " << diff << "MiB ";
    }
  }

 private:
  std::string m_name;
  size_t m_origRss;
};

//! stl allocator using malloc/free
template <typename T>
class malloc_allocator : public std::allocator<T> {
 public:
  typedef T* pointer;
  typedef size_t size_type;

  template <typename U>
  struct rebind {
    typedef malloc_allocator<U> other;
  };

  pointer allocate(size_type n, const void* hint = 0) {
    if (!n) return NULL;

    pointer res = static_cast<pointer>(std::malloc(n * sizeof(T)));
    if (!res) throw std::bad_alloc();
    return res;
  }

  void deallocate(pointer p, size_type n) { std::free(p); }

  malloc_allocator() throw() : std::allocator<T>() {}

  template <typename U>
  malloc_allocator(const malloc_allocator<U>& a) throw()
      : std::allocator<T>(a) {}

  ~malloc_allocator() throw() {}
};

//! Keeps track of allocations. Does nothing while OVERLOAD_ALLOCATORS is not
//! set in .cc file
struct AllocStats {
  typedef const void* alloc_t;

  static void print(std::ostream& out, bool sortByCalls = false,
                    bool withDetails = true, bool thenReset = false);

  static void ack_malloc(alloc_t data, size_t size);
  static void ack_free(alloc_t data);

  //! Whether to count only non-freed memory or total allocations
  static bool s_count_total_allocated;

 private:
  static void print_size(size_t size, std::ostream& out);

  typedef std::pair<void*, void*> caller_address;
  typedef std::pair<size_t, unsigned> caller_stats;

  typedef std::map<
      caller_address, caller_stats, std::less<caller_address>,
      malloc_allocator<std::pair<const caller_address, caller_stats> > >
      CallersType;

  typedef std::map<
      alloc_t, std::pair<caller_address, size_t>, std::less<alloc_t>,
      malloc_allocator<
          std::pair<const alloc_t, std::pair<caller_address, size_t> > > >
      AllocationsType;

  static MutexType s_mutex;
  static AllocationsType s_allocations;
  static CallersType s_callers;
};

struct compare_c_strings
    : public std::binary_function<const char*, const char*, bool> {
  bool operator()(const char* s1, const char* s2) const {
    return std::strcmp(s1, s2) < 0;
  }
};

//! Utility class using /proc/self/maps and /bin/addr2line to translate
//! addresses to source file lines
class Addr2Line {
 public:
  Addr2Line();
  const char* get_file_and_line(size_t addr, const char* libname);

  static char* demangle(char* mangled_name);

 private:
  typedef std::map<const char*, size_t, compare_c_strings,
                   malloc_allocator<std::pair<const char* const, size_t> > >
      AddressesMap;

  static AddressesMap s_addresses;

  static std::string get_library_base_name(const char* libname);
  static size_t get_library_address(const AddressesMap& addresses,
                                    const char* libname);
  static void get_library_addresses(AddressesMap& addresses,
                                    bool verbose = false);

  static const size_t npos;
};

}  // namespace strandsim

#endif  // STRANDSIM_MEMORY_HH
