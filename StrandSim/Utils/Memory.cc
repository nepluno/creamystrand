/**
 * \copyright 2012 Gilles Daviet
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Memory.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "ThreadUtils.hh"
#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <cxxabi.h>
#include <execinfo.h>
#endif
#include <fcntl.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

//#define OVERLOAD_ALLOCATORS

namespace strandsim {

Addr2Line::AddressesMap Addr2Line::s_addresses;
const size_t Addr2Line::npos = (size_t)-1;

bool AllocStats::s_count_total_allocated = false;
AllocStats::AllocationsType AllocStats::s_allocations;
AllocStats::CallersType AllocStats::s_callers;
MutexType AllocStats::s_mutex;

void AllocStats::ack_malloc(alloc_t data, size_t size) {
  static const int target_frame = 4;
  void* addrlist[target_frame];

#ifndef WIN32
  if (backtrace(addrlist, target_frame) < target_frame) {
    return;
  }
#endif

  caller_address caller =
      std::make_pair(addrlist[target_frame - 2], addrlist[target_frame - 1]);
  ;

  LockGuard lock(s_mutex);
  if (!s_count_total_allocated) {
    s_allocations[data] = std::make_pair(caller, size);
  }
  caller_stats& stats = s_callers[caller];
  stats.first += size;
  ++stats.second;
}

void AllocStats::ack_free(alloc_t data) {
  if (s_count_total_allocated) return;

  LockGuard lock(s_mutex);

  auto alloc_it = s_allocations.find(data);
  if (alloc_it != s_allocations.end()) {
    s_callers[alloc_it->second.first].first -= alloc_it->second.second;
    s_allocations.erase(alloc_it);
  }
}

struct memory_record {
  struct stack {
    void* addr;
    const char* lib_name;
    size_t size;
    unsigned calls;
  };

  const char* func_name;
  size_t size;
  size_t order;
  std::list<stack, malloc_allocator<stack> > stacks;

  bool operator<(const memory_record& other) const {
    return order > other.order;
  }
};

void AllocStats::print_size(size_t size, std::ostream& out) {
  if (size >> 32) {
    out << (size >> 30) << " GiB";
  } else if (size >> 22) {
    out << (size >> 20) << " MiB";
  } else if (size >> 12) {
    out << (size >> 10) << " kiB";
  } else {
    out << (size) << " bytes";
  }
}

void AllocStats::print(std::ostream& out, bool sortByCalls, bool withDetails,
                       bool thenReset) {
#ifndef WIN32
#ifndef OVERLOAD_ALLOCATORS
  out << " AllocStats::print requires OVERLOAD_ALLOCATORS to be defined in "
         "Memory.cc\n";
  return;
#endif

  std::list<memory_record, malloc_allocator<memory_record> > records;
  std::vector<char*, malloc_allocator<char*> > func_names;
  char **symbollist = NULL, **symbollist_backup = NULL;

  {
    LockGuard lock(s_mutex);

    const unsigned nAllocs = s_callers.size();
    func_names.resize(nAllocs);

    std::vector<void*, malloc_allocator<void*> > symbols(nAllocs);
    std::vector<void*, malloc_allocator<void*> > symbols_backup(nAllocs);
    std::vector<caller_stats, malloc_allocator<caller_stats> > stats(nAllocs);

    std::map<char*, memory_record, compare_c_strings,
             malloc_allocator<std::pair<char* const, memory_record> > >
        symbol_names;

    unsigned i = 0;
    for (auto alloc_it = s_callers.begin(); alloc_it != s_callers.end();
         ++alloc_it) {
      symbols[i] = alloc_it->first.first;
      symbols_backup[i] = alloc_it->first.second;
      stats[i] = alloc_it->second;
      ++i;
    }
    symbollist = backtrace_symbols(&symbols[0], nAllocs);
    symbollist_backup = backtrace_symbols(&symbols_backup[0], nAllocs);

    for (unsigned i = 0; i < nAllocs; ++i) {
      if ((func_names[i] =
               Addr2Line::demangle(symbollist[i])))  // <- assignment
      {
        memory_record::stack stack;
        stack.size = stats[i].first;
        stack.calls = stats[i].second;
        stack.addr = symbols[i];
        stack.lib_name = symbollist[i];

        if (strncmp(func_names[i], "std::", 5) == 0 ||
            strncmp(func_names[i], "Eigen::", 7) == 0) {
          // Standard lib function, not that useful
          free(func_names[i]);
          func_names[i] = Addr2Line::demangle(symbollist_backup[i]);

          stack.addr = symbols_backup[i];
          stack.lib_name = symbollist_backup[i];
        }

        memory_record& rec = symbol_names[func_names[i]];
        rec.size += stack.size;
        if (sortByCalls) {
          rec.order += stack.calls;
        } else {
          rec.order += stack.size;
        }
        rec.func_name = func_names[i];
        rec.stacks.push_back(stack);
      }
      //        symbol_names[ symbollist[i] ] = sizes[i] ;
    }

    for (auto symb_it = symbol_names.begin(); symb_it != symbol_names.end();
         ++symb_it) {
      if (symb_it->second.size) {
        records.push_back(symb_it->second);
      }
    }

    if (thenReset) {
      s_allocations.clear();
      s_callers.clear();
    }
  }
  records.sort();

  unsigned max = 16;
  unsigned numPrinted = 0;
  for (auto symb_it = records.begin();
       symb_it != records.end() && numPrinted < max; ++symb_it) {
    if (!symb_it->size) break;

    out << "[" << (++numPrinted) << "] ";
    print_size(symb_it->size, out);
    out << " have been allocated by \n" << symb_it->func_name << std::endl;
    if (withDetails) {
      for (auto stack_it = symb_it->stacks.begin();
           stack_it != symb_it->stacks.end(); ++stack_it) {
        if (stack_it->size || (sortByCalls && stack_it->calls)) {
          out << "  +-> ";
          print_size(stack_it->size, out);
          out << " in " << stack_it->calls << " calls by "
              << Addr2Line().get_file_and_line((size_t)stack_it->addr,
                                               stack_it->lib_name);
        }
      }
    }
  }

  for (unsigned i = 0; i < func_names.size(); ++i) {
    if (func_names[i]) {
      free(func_names[i]);
    }
  }
  if (symbollist) free(symbollist);
  if (symbollist_backup) free(symbollist_backup);
#endif
}

size_t getResidentMemory() {
  size_t rss = 0;
#ifdef __APPLE__
  std::cerr << "/proc not supported on OS X. " << std::endl;
#else

  // 'file' stat seems to give the most reliable results
  //
  std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);

  // dummy vars for leading entries in stat that we don't care about
  //
  std::string pid, comm, state, ppid, pgrp, session, tty_nr;
  std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  std::string utime, stime, cutime, cstime, priority, nice;
  std::string O, itrealvalue, starttime, vsize;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
      tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
      stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
      starttime >> vsize >> rss;  // don't care about the rest

  stat_stream.close();
#ifdef WIN32
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  rss = (rss * si.dwPageSize) >> 20;  // In MB
#else
  rss = (rss * sysconf(_SC_PAGE_SIZE)) >> 20;  // In MB
#endif

#endif

  return rss;
}

// Retrieves the file name from a path
std::string Addr2Line::get_library_base_name(const char* libname) {
  const char* last_slash = std::strrchr(libname, '/');
  return std::string(last_slash ? last_slash + 1 : libname);
}

// Retrieves previously found address from a library path
size_t Addr2Line::get_library_address(const AddressesMap& addresses,
                                      const char* libname) {
  const auto& it = addresses.find(libname);
  if (it != addresses.end()) {
    return it->second;
  }

  return npos;
}

// Reads all dynamic libraries base addresses from /proc/PID/smaps
void Addr2Line::get_library_addresses(AddressesMap& addresses, bool verbose) {
  //    char path[256];
  //    snprintf(path, sizeof path, "/proc/%d/smaps", getpid());
  const char* path = "/proc/self/maps";

  if (verbose)
    std::cerr << "Reading dynamic libraries addresses from " << path << "... "
              << std::endl;

  const int addr_max = 17;  // Including terminal \0
  char addr[addr_max];
  char lib_name[PATH_MAX];

  const int buf_size = 96;
  char buf[buf_size];

  FILE* maps = fopen(path, "r");
  if (maps) {
    bool line_done = true;
    const char* last_space = NULL;
    char *cur_buf, *cur_libname_pos;

    while (!feof(maps) && fgets(buf, buf_size, maps)) {
      cur_buf = buf;

      if (line_done) {
        cur_libname_pos = lib_name;

        // Beginning a new line, read address

        cur_buf = strchr(buf, '-');
        if (!cur_buf) continue;

        const int addr_len = cur_buf - buf;
        if (addr_len >= addr_max) continue;

        strncpy(addr, buf, addr_len + 1);
        addr[addr_len] = '\0';

        *cur_buf++ = '\0';
      }
      line_done = false;

      last_space = cur_buf;
      while (*cur_buf) {
        if (*cur_buf == ' ' || *cur_buf == '\t') {
          cur_libname_pos = lib_name;
          last_space = cur_buf + 1;
        } else if (*cur_buf == '\n') {
          line_done = true;
          break;
        }
        ++cur_buf;
      }

      const int n_since_space = cur_buf - last_space;
      if (n_since_space) {
        // Append to lib_name
        strncpy(cur_libname_pos, last_space, n_since_space);
        cur_libname_pos += n_since_space;
      }

      // Line finished, save address
      if (line_done) {
        *cur_libname_pos = '\0';
        if (addresses.find(lib_name) == addresses.end()) {
          if (verbose)
            std::cerr << " " << lib_name << " at " << addr << std::endl;

          // Allocate new string that will never be freed, but we dont' care
          // as this function will only be called once
          size_t len = strlen(lib_name);
          char* lib_name_alloc = (char*)malloc(len + 1);
          strncpy(lib_name_alloc, lib_name, len + 1);

          addresses[lib_name_alloc] = std::strtoul(addr, NULL, 16);
        }
      }
    }

    fclose(maps);
  }

  if (verbose) std::cerr << "done. " << std::endl;
}

// Uses addr2line to find the source line for a given memory address
const char* Addr2Line::get_file_and_line(size_t addr, const char* libname) {
  // storage array for absolute library path
  static char s_absolute_path[PATH_MAX];
  // storage array for file name + line number
  static const unsigned s_filename_max_size = 256;
  static char s_filename[s_filename_max_size];

  const char* ret = NULL;
#ifdef WIN32
  if (GetFullPathName(libname, sizeof(s_absolute_path), s_absolute_path, NULL))
#else
  if (realpath(libname, s_absolute_path))
#endif
  {
    libname = s_absolute_path;
  }

  const unsigned libname_length = strlen(libname);

  size_t base_addr = get_library_address(s_addresses, libname);
  if (base_addr == npos) return libname;

  // Assembles addr2line command
  char* cmd = (char*)std::malloc(libname_length + 128);
  sprintf(cmd, "/usr/bin/addr2line -s -e \"%s\" %#lx", libname,
          addr - base_addr);
//    fprintf( stderr, "[ %s ]\n", cmd ) ;
#ifdef WIN32
  FILE* f = _popen(cmd, "r");
#else
  FILE* f = popen(cmd, "r");
#endif

  if (f) {
    if (!feof(f) && fgets(s_filename, s_filename_max_size, f)) {
      // 5 chars or less means that addr2line failed
      if (strlen(s_filename) > 5) {
        ret = s_filename;
      }
    }
#ifdef WIN32
    _pclose(f);
#else
    pclose(f);
#endif
  }

  free(cmd);

  return ret;
}

char* Addr2Line::demangle(char* mangled_name) {
  char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

  // find parentheses and +address offset surrounding the mangled name:
  // ./module(function+0x15c) [0x8048a6d]
  for (char* p = mangled_name; *p; ++p) {
    if (*p == '(')
      begin_name = p;
    else if (*p == '+')
      begin_offset = p;
    else if (*p == ')' && begin_offset) {
      end_offset = p;
      break;
    }
  }

  if (begin_name && begin_offset && end_offset && begin_name < begin_offset) {
    *begin_offset++ = '\0';
    *begin_name++ = '\0';

    int status;
#ifdef WIN32
    return begin_name;
#else
    return abi::__cxa_demangle(begin_name, NULL, NULL, &status);
#endif
  }

  return NULL;
}

Addr2Line::Addr2Line() {
  if (s_addresses.empty()) {
    get_library_addresses(s_addresses);
  }
}

}  // namespace strandsim

#ifdef OVERLOAD_ALLOCATORS

void* operator new(size_t size) throw(std::bad_alloc) {
  void* ret = std::malloc(size);
  if (!ret) throw std::bad_alloc();
  strandsim::AllocStats::ack_malloc(ret, size);
  return ret;
}

void* operator new[](size_t size) throw(std::bad_alloc) {
  void* ret = std::malloc(size);
  if (!ret) throw std::bad_alloc();
  strandsim::AllocStats::ack_malloc(ret, size);
  return ret;
}

void operator delete(void* data) throw() {
  std::free(data);
  strandsim::AllocStats::ack_free(data);
}

void operator delete[](void* data) throw() {
  std::free(data);
  strandsim::AllocStats::ack_free(data);
}

#endif
