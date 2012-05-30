#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <boost/crc.hpp>
#include <MyMath.h>
#include <vector>


// Please note the MD5 Copyright notice in the 
// implementation file IO.cc!

namespace io {

  class CheckedIo;

  std::ostream& operator <<(std::ostream& os, CheckedIo &io);

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Helper class to write complex numbers to a std::ostream and
  ///  calculate the md5 checksum on-the-fly.
  ///
  ///  \author Dirk Hesse <herr.dirk.hesse@gmail.com>
  ///  \date Fri May 25 16:25:45 2012

  class CheckedIo {
  public:
    CheckedIo() : os("DATA_OUT", std::ios::ate | std::ios::app |
                     std::ios::binary), bcount(0), buffcnt(0) { 
      h[0] = 0x67452301;
      h[1] = 0xEFCDAB89;
      h[2] = 0x98BADCFE;
      h[3] = 0x10325476;
    } 
    ~CheckedIo() {
      os.close();
      finalize();
      std::cout << "md5 checksum: " << *this << std::endl;
    }
    void finalize();
    void write(const Cplx &c){
      //w[buffcnt] = c.re;
      //w[buffcnt+1].d = c.im;
      reinterpret_cast<double *>(w)[buffcnt] = c.re;
      reinterpret_cast<double *>(w)[buffcnt+1] = c.im;
      buffcnt += 2;
      if (buffcnt == 8){
        md5process();
        buffcnt = 0;
        ++bcount;
        os.write(reinterpret_cast<char const*>(w), 16*sizeof(unsigned));
      }
    }
    std::vector<unsigned> get_h() const {
      return std::vector<unsigned>(h, h+4);
    }
  private:
    std::ofstream os;
    unsigned h[4];
    static const unsigned r[], k[];
    unsigned w[16];
    // block count
    unsigned bcount;
    // buffer count
    int buffcnt;
    void md5process();
    unsigned lrol(const unsigned& u, const unsigned & shift) {
      return (u << shift) | (u >> (32 - shift));
    }
  };
}


#endif /* _IO_H_ */
