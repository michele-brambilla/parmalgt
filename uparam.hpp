#ifndef UPARAM_H
#define UPARAM_H

#include <string>
#include <algorithm>
#include <map>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///
///  Minimalist framework to read/write parameters.
///
///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
///  \date Wed May 30 16:39:24 2012

namespace uparam {
  class Param {
  public:
    explicit Param (const std::string& filename) : fname(filename) { }
    void read(){
      std::ifstream inf(fname.c_str());
      if (!inf.is_open()) {
        std::cerr << "PARAMETER FILE " << fname  << " NOT FOUND!\n";
        throw std::exception();
      }
      while (inf.good()) {
        std::string name = get_str(inf);
        params[name] = get_str(inf);
      }
      inf.close();
    }
    void write(){
      std::map<std::string, std::string>::const_iterator i;
      std::ofstream of(fname.c_str(), std::ios::trunc);
      for (i = params.begin(); i != params.end(); ++i)
        of << i->first << "  " << i->second << "\n";
      of.close();
    }
    std::string& operator[](const std::string& s){
      return params[s];
    }
  private:
    std::string fname;
    std::map<std::string, std::string> params;
    std::string to_upper(std::string& in){
      std::transform(in.begin(), in.end(), in.begin(), 
                     (int(*)(int))std::toupper);
      return in;
    }
    std::string get_str(std::ifstream& in) {
      std::string next;
      in >> next;
      char dummy[256];
      while (next[0] == '#' && in.good()) {
        in.getline(dummy, 256);
        in >> next;
      }
      if (next[0] == '#')
        next = "";
      return next;
    }
  };
}
#endif
