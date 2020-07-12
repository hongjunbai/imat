#ifndef SPLIT_STR_H
#define SPLIT_STR_H
#include <string>
#include <vector>
#include "macros.h"

/* USAGE:
 * Given:
 * string a("Error! splite_str() is not defined for T != std::string!");
 * std::vector<std::string> list;
 * eg. 1
 * split_str(&list, a); // default set space as " \t\n\r"
 * eg. 2
 * std::string b=" ";
 * split_str(&list, a, b);
 * eg. 3
 * split_str(&list, a, string("# "));
 *
 * BEWARE! The following code does not work.
 * split_str(&list, a, "# ");
 */

// It seems that inline helps to avoid multiple definition while compiling
template <class T>
inline void split_str(std::vector<T> *list,
               const T &source,
               const T &delims = " \t\n\r") {
  MESHOW("Error! splite_str() is not defined for T != std::string!");
}

template <>
inline void split_str(std::vector<std::string> *list,
               const std::string &source,
               const std::string &delims) {
  std::string str;
  list->clear();
  std::string::size_type pos=0, len=source.size();
  while (pos<len){
    str = "";
    // Remove any delimiters 
    while ((delims.find(source[pos]) != std::string::npos) && (pos<len)){
      ++pos;
    }
    // Save token data
    while ((delims.find(source[pos]) == std::string::npos) && (pos<len)){
      str += source[pos++];
    }
    // Put valid str buffer into the supplied list
    if(! str.empty()) list->push_back(str);
  }
}

template <class T>
inline void strip_str(T *source, const T &delims = " \t\n\r") {
  MESHOW("Error! strip_str() is not defined for T != std::string!");
}

template <>
inline void strip_str(std::string *source, const std::string &delims) {
  std::string::size_type len=source->size();
  std::string::size_type cnt0=0, idx=len-1;
  while ((delims.find((*source)[cnt0]) != std::string::npos) && (cnt0<len)){
    ++cnt0;
  }
  while ((delims.find((*source)[idx]) != std::string::npos) && (idx>=0)){
    --idx;
  }
  if (cnt0 > idx) {
    source->clear();
  } else {
    source->erase(source->begin(), source->begin() + cnt0);
    source->erase(source->end()-(len-1 - idx), source->end());
  }
}

#endif // SPLIT_STR_H
