// -*- C++ -*-
//    Title: common.h
//  Created: Mon Dec 17 20:25:17 2018
//   Author: Tim Weyrich <t.weyrich@cs.ucl.ac.uk>
//      $Id: $
//
// copyright (c) 2017, University College London
//

#ifndef __COMMON_H__
#define __COMMON_H__

#include <assert.h>

#include <string>
#include <iterator>

//////////////////////////////////////////////////////////////////////////////
// Helpers
//

template<typename T>
std::string  str(const T &val) {
    return std::to_string(val);
}

template<typename T>
std::string  str(const glm::tvec3<T> &val) {
    char  s[1024];
    sprintf(s, "[%g; %g; %g]", val.x, val.y, val.z);  // output follows Matlab convention for easier copy-and-paste
    return std::string(s);
}

template<typename V>
std::string  str(const std::vector<V> &v) {
    std::string  s;
    s.append("[ ");
    for (int i=0; i<(int)v.size(); ++i) {
        if (i > 0)
            s.append(", ");
        s.append(str(v[i]));
    }
    s.append(" ]");
    return s;
}

template<typename T>
bool  anyIsnan(const T &val) {
    return std::isnan(val);
}

template<typename T>
bool  anyIsnan(const glm::tvec3<T> &v) {
    return (std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z));
}

template<typename V>
bool  anyIsnan(const std::vector<V> &v) {
    for (int i=0; i<(int)v.size(); ++i)
        if (anyIsnan(v[i]))
            return true;
    return false;
}

#define  assertNoNan(x)  do{assertNoNanFun((x),__FILE__,__LINE__,#x);}while(0)

template<typename T>
void  assertNoNanFun(const T &val, const char *file, int line, const char *expr) {
    if (anyIsnan(val)) {
        std::cerr << file << ":" << line << ": NaN contained in `" << expr << "':\n" << str(val) << std::endl;
        assert(!"isnan test failed"); // let's break in a way that allow for debugging from here
    }
}

inline double  randSigned() { return 2.0*rand()/RAND_MAX-1.0; }

inline glm::dvec3  randDirection3()
{
    glm::dvec3  dir;
    double  len;
    do {
         dir = glm::dvec3(randSigned(), randSigned(), randSigned());
         len = glm::length(dir);
    } while (len > 1);
    return (1.0/len) * dir;
}

inline glm::dmat3  randRotationMatrix3()
{
    auto  u = randDirection3();
    auto  v = glm::normalize(glm::cross(randDirection3(), u));
    auto  w = glm::cross(u, v);
    return glm::dmat3(u,v,w);
}

template<typename C>
bool  anyRepeat(const C &container) {
    for (auto it=container.begin(); it!=container.end(); ++it)
        if (glm::length(*it - *std::next(it)) <= 1e-6)
            return true;
    return container.size() < 2 || glm::length(container.front() - container.back()) <= 1e-6;
}

#define  assertNoRepeat(x)  do{assertNoRepeatFun((x),__FILE__,__LINE__,#x);}while(0)

template<typename T>
void  assertNoRepeatFun(const T &container, const char *file, int line, const char *expr) {
    if (anyRepeat(container)) {
        std::cerr << file << ":" << line << ": repeated value in `" << expr << "':\n" << str(container) << std::endl;
        assert(!"no-repeat test failed"); // let's break in a way that allow for debugging from here
    }
}

#endif
