#include "IO.hpp"

// The CheckedIo class uses the RSA Data Security, Inc. MD5
// Message-Digest Algorithm.
//
// Copyright (C) 1991-2, RSA Data Security, Inc. Created 1991. All
// rights reserved.
//
// License to copy and use this software is granted provided that it
// is identified as the "RSA Data Security, Inc. MD5 Message-Digest
// Algorithm" in all material mentioning or referencing this software
// or this function.
//
// License is also granted to make and use derivative works provided
// that such works are identified as "derived from the RSA Data
// Security, Inc. MD5 Message-Digest Algorithm" in all material
// mentioning or referencing the derived work.
//
// RSA Data Security, Inc. makes no representations concerning either
// the merchantability of this software or the suitability of this
// software for any particular purpose. It is provided "as is"
// without express or implied warranty of any kind.
//
// These notices must be retained in any copies of any part of this
// documentation and/or software.


/*
 * static members
 */
const unsigned io::CheckedIo::r[] = { 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17,
    22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
    4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15, 21,
    6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21 };
const unsigned io::CheckedIo::k[] = { 3614090360u, 3905402710u, 606105819u,
    3250441966u, 4118548399u, 1200080426u, 2821735955u, 4249261313u,
    1770035416u, 2336552879u, 4294925233u, 2304563134u, 1804603682u,
    4254626195u, 2792965006u, 1236535329u, 4129170786u, 3225465664u,
    643717713u, 3921069994u, 3593408605u, 38016083u, 3634488961u, 3889429448u,
    568446438u, 3275163606u, 4107603335u, 1163531501u, 2850285829u,
    4243563512u, 1735328473u, 2368359562u, 4294588738u, 2272392833u,
    1839030562u, 4259657740u, 2763975236u, 1272893353u, 4139469664u,
    3200236656u, 681279174u, 3936430074u, 3572445317u, 76029189u, 3654602809u,
    3873151461u, 530742520u, 3299628645u, 4096336452u, 1126891415u,
    2878612391u, 4237533241u, 1700485571u, 2399980690u, 4293915773u,
    2240044497u, 1873313359u, 4264355552u, 2734768916u, 1309151649u,
    4149444226u, 3174756917u, 718787259u, 3951481745u };

void io::CheckedIo::finalize() {

  // get the length of the last io
  unsigned last_io = buffcnt*sizeof(double);
  /*
   * the algorithm ends as follows:
   * o append the bit "1" to the message
   * o pad the message with "0", until length(message)[bit] % 512 == 448
   * o append length(message)[bit] (before padding) to message (as 64 bit unsigned)
   */
  unsigned fsize[] = { 0, 0 };
  // handle files of more than 4GB
  const unsigned mod = 536870912; // max_unsigned_value / 8
  fsize[0] = (bcount * 64 + last_io) % mod;
  fsize[1] = (bcount * 64 + last_io - fsize[0]) / mod;
  for (int i = 0; i < 2; i++)
    fsize[i] *= 8;
  // pointer in w (unsigned char) to the last read position:
  int ptr = last_io;
  // do the padding
  w.c[ptr++] = 128;
  while (ptr % 64 != 56) {
    if (ptr >= 64) {
      md5process();
      ptr = 0;
    }
    w.c[ptr++] = 0;
  }
  w.u[14] = fsize[0];
  w.u[15] = fsize[1];
  md5process();

}

void io::CheckedIo::md5process() {
  // the results so far
  unsigned a = h[0], b = h[1], c = h[2], d = h[3];
  unsigned f, g, tmp;
  for (int i = 0; i < 16; i++) {
    f = (b & c) | ((~b) & d);
    g = i;
    tmp = d;
    d = c;
    c = b;
    b += lrol(a + f + k[i] + w.u[g], r[i]);
    a = tmp;
  }
  for (int i = 16; i < 32; i++) {
    f = (d & b) | ((~d) & c);
    g = (5 * i + 1) % 16;
    tmp = d;
    d = c;
    c = b;
    b += lrol(a + f + k[i] + w.u[g], r[i]);
    a = tmp;
  }
  for (int i = 32; i < 48; i++) {
    f = b ^ c ^ d;
    g = (3 * i + 5) % 16;
    tmp = d;
    d = c;
    c = b;
    b += lrol(a + f + k[i] + w.u[g], r[i]);
    a = tmp;
  }
  for (int i = 48; i < 64; i++) {
    f = c ^ (b | (~d));
    g = (7 * i) % 16;
    tmp = d;
    d = c;
    c = b;
    b += lrol(a + f + k[i] + w.u[g], r[i]);
    a = tmp;
  }
  h[0] += a;
  h[1] += b;
  h[2] += c;
  h[3] += d;
}

std::ostream& io::operator <<(std::ostream& os, CheckedIo &io) {
    //unsigned char* w = reinterpret_cast<unsigned char*> (md5.h);
    std::vector<unsigned> h = io.get_h();
    unsigned char* w = reinterpret_cast<unsigned char*> (&(h[0]));
    os << std::hex;
    for (int i = 0; i < 16; i++)
      os << (int) w[i];
    os << std::dec;
    return os;
}
