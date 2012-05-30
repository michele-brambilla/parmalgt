#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <MyMath.h>
#include <vector>
#include <sstream>
#include <uparam.hpp>


// Note the MD5 Copyright notice in the 
// implementation file IO.cc!

namespace io {

  // The data on which md5 operates has to be accessed as double,
  // unsigned and unsigned char. This union will help here...
  
  namespace detail {
    // NOTE THAT THESE ARE ARCHITECTURE DEPENDENT!
    const int unsigned_per_double = 2;
    const int char_per_double = 8;

    union md5atom {
      double d[8];
      unsigned u[8 * unsigned_per_double];
      unsigned char c[8 * char_per_double];
    };
  }

  class CheckedIo;

  std::ostream& operator <<(std::ostream& os, CheckedIo &io);

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Helper calculate the md5 checksum on-the-fly.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Fri May 25 16:25:45 2012

  class CheckedIo {
  public:
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Constructor.
    ///
    ///  Initialize the md5 algorithm.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed May 30 18:37:52 2012
    CheckedIo() : bcount(0), buffcnt(0) { 
      h[0] = 0x67452301;
      h[1] = 0xEFCDAB89;
      h[2] = 0x98BADCFE;
      h[3] = 0x10325476;
    } 
    ~CheckedIo() { }
    void finalize();
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Return the md5 as a vector.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed May 30 18:38:30 2012
    std::vector<unsigned> get_h() const {
      return std::vector<unsigned>(h, h+4);
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Return the md5 as a string.
    ///
    ///  Note: This calls finialize, so it should only be called ONCE!
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed May 30 18:38:39 2012
    std::string md5(){
      finalize();
      unsigned char* w = reinterpret_cast<unsigned char*> (h);
      std::stringstream s;
      s << std::hex;
      for (int i = 0; i < 16; i++)
        s << (int) w[i];
      s << std::dec;
      return s.str();
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Process a dobule.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed May 30 18:40:06 2012
    void process(const double& d){
      w.d[buffcnt++] = d;
      if (buffcnt == 8){
        md5process();
        buffcnt = 0;
        ++bcount;
      }
    }
  private:
    unsigned h[4];
    static const unsigned r[], k[];
    detail::md5atom w;
    // block count
    unsigned bcount;
    // buffer count
    int buffcnt;
    void md5process();
    unsigned lrol(const unsigned& u, const unsigned & shift) {
      return (u << shift) | (u >> (32 - shift));
    }
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Writes complex numbers to a file and calculates the md5
  ///  checksum on-the-fly.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed May 30 18:40:15 2012
  class CheckedOut : public CheckedIo {
  public:
    CheckedOut() : CheckedIo(), os("gauge.cfg", std::ios::trunc |
                                   std::ios::binary ),
                   param("gauge.info"){ }
    ~CheckedOut() { 
      os.close();
      param["md5"] = CheckedIo::md5();
      param.write();
    }
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    ///
    ///  Write a complex to disk.
    ///
    ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
    ///  \date Wed May 30 18:40:45 2012
    void write(const Cplx &c){
      CheckedIo::process(c.re);
      CheckedIo::process(c.im);
      os.write(reinterpret_cast<char const*>(&(c.re)), sizeof(double));
      os.write(reinterpret_cast<char const*>(&(c.im)), sizeof(double));
    }
  private:
    std::ofstream os;
    uparam::Param param;
  };

  class IoError : public std::exception {};


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Read complex numbers from disk.
  ///
  ///  Check if the md5 checksum is right.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Wed May 30 18:40:58 2012
  class CheckedIn : public CheckedIo {
  public:
    CheckedIn() : CheckedIo(), is("gauge.cfg", std::ios::binary ),
                  param("gauge.info"){
      param.read();
    }
    ~CheckedIn() {
      is.close();
      std::string md5 = CheckedIo::md5();
      if (md5 != param["md5"]){
        std::cout << "'" << md5 << "'" << std::endl;
        std::cout << "'" << param["md5"] << "'" << std::endl;
        throw IoError();
      }
    }
    void read(Cplx& c){
      is.read(reinterpret_cast<char*>(&(c.re)), sizeof(double));
      is.read(reinterpret_cast<char*>(&(c.im)), sizeof(double));
      CheckedIo::process(c.re);
      CheckedIo::process(c.im);
    }
  private:
    std::ifstream is;
    uparam::Param param;
  };
}


#endif /* _IO_H_ */
