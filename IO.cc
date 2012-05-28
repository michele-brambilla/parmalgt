#include <IO.hpp>
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

  // count of the blocks
  //unsigned bcount = 0;

  //while (!f->fail()) {
  //  // get 512 bit chunk from file
  //  f->read(reinterpret_cast<char*> (w), 64);
  //  if (!f->fail()) {
  //    md5process();
  //    bcount++;
  //  }
  //}
  // get the length of the last read
  unsigned last_read = buffcnt*sizeof(double);
  /*
   * the algorithm ends as follows:
   * o append the bit "1" to the message
   * o pad the message with "0", until length(message)[bit] % 512 == 448
   * o append length(message)[bit] (before padding) to message (as 64 bit unsigned)
   */
  unsigned fsize[] = { 0, 0 };
  // handle files of more than 4GB
  const unsigned mod = 536870912; // max_unsigned_value / 8
  fsize[0] = (bcount * 64 + last_read) % mod;
  fsize[1] = (bcount * 64 + last_read - fsize[0]) / mod;
  for (int i = 0; i < 2; i++)
    fsize[i] *= 8;
  // get a finer resolution on w
  unsigned char* wc = reinterpret_cast<unsigned char*> (w);
  // pointer in wc to the last read position:
  int ptr = last_read;
  std::cout << last_read << "\n";
  // do the padding
  wc[ptr++] = 128;
  while (ptr % 64 != 56) {
    if (ptr >= 64) {
      md5process();
      ptr = 0;
    }
    wc[ptr++] = 0;
  }
  //w[14] = fsize[0];
  //w[15] = fsize[1];
  w[7].u[0] = fsize[0];
  w[7].u[1] = fsize[1];
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
    b += lrol(a + f + k[i] + w[g/2].u[g%2], r[i]);
    a = tmp;
  }
  for (int i = 16; i < 32; i++) {
    f = (d & b) | ((~d) & c);
    g = (5 * i + 1) % 16;
    tmp = d;
    d = c;
    c = b;
    b += lrol(a + f + k[i] + w[g/2].u[g%2], r[i]);
    a = tmp;
  }
  for (int i = 32; i < 48; i++) {
    f = b ^ c ^ d;
    g = (3 * i + 5) % 16;
    tmp = d;
    d = c;
    c = b;
    b += lrol(a + f + k[i] + w[g/2].u[g%2], r[i]);
    a = tmp;
  }
  for (int i = 48; i < 64; i++) {
    f = c ^ (b | (~d));
    g = (7 * i) % 16;
    tmp = d;
    d = c;
    c = b;
    b += lrol(a + f + k[i] + w[g/2].u[g%2], r[i]);
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
    return os;
}
