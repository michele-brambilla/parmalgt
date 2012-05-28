#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <boost/crc.hpp>
#include <MyMath.h>
#include <vector>


namespace io {

  namespace detail {
    union md5atom {
      double d;
      unsigned u[sizeof(double)]; // should be unsigned[2] on x86...
    };
  }

  class CheckedIo;

  std::ostream& operator <<(std::ostream& os, CheckedIo &io);

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///
  ///  Helper class to write complex numbers to a std::ostream and
  ///  calculate the crc32 checksum on-the-fly.
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
      //std::cout << "Checksum: " << crc.checksum() << std::endl;
    }
    void finalize();
    void write(const Cplx &c){
      w[buffcnt].d = c.re;
      w[buffcnt+1].d = c.im;
      buffcnt += 2;
      if (buffcnt == 8){
        md5process();
        buffcnt = 0;
        ++bcount;
      }
      //os.write(reinterpret_cast<char const*>(&(c.re)), sizeof(double));
      //os.write(reinterpret_cast<char const*>(&(c.im)), sizeof(double));
      //crc.process_block(reinterpret_cast<void const*>((&c.re)),
      //                  reinterpret_cast<void const*>((&c.re + 1)));
      //crc.process_block(reinterpret_cast<void const*>((&c.im)),
      //                  reinterpret_cast<void const*>((&c.im + 1)));
    }
    std::vector<unsigned> get_h() const {
      return std::vector<unsigned>(h, h+4);
    }
  private:
    std::ofstream os;
    //boost::crc_optimal<32> crc;
    unsigned h[4];
    static const unsigned r[], k[];
    detail::md5atom w[8];
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
