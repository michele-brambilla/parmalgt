#ifndef _COMMUNICATOR_H
#define _COMMUNICATOR_H

#include <Background.h>
#include <Geometry.hpp>
#include <Point.hpp>
#include <newQCDpt.h>
#include <IO.hpp>

#include <mpi.h>

#include <map>
#include <iostream>
#include <fstream>
#include <uparam.hpp>
#include <stdlib.h>



namespace kernels {
  template <class Field_t>
  struct base_types;
};


namespace comm {

  template<class Field_t>
  struct Communicator : public kernels::base_types<Field_t> {

    typedef typename kernels::base_types<Field_t>::data_t data_t;
    typedef typename kernels::base_types<Field_t>::point_t Point;
    typedef typename kernels::base_types<Field_t>::direction_t Direction;
    static const int ORD = kernels::base_types<Field_t>::order;
    static const int DIM = kernels::base_types<Field_t>::n_dim;

    typedef typename array_t<std::pair<int,int>, DIM>::Type nb_t;
    typedef typename geometry::Geometry<DIM> Geometry;
    typedef typename Geometry::extents_t extents_t;
    typedef typename std::vector< std::pair<std::vector<data_t>,
					    std::vector<data_t> > > buffer_t;
    typedef typename std::vector<data_t>::iterator iterator_t;

    /// Buffer for communication
    buffer_t send_buff_;
    buffer_t recv_buff_;

    static nb_t nb_;
    static int numprocs_;
    static int rank_;

    Communicator(const extents_t& e) {
      for (int i = 0; i < DIM; ++i){
	send_buff_.push_back
	  (std::make_pair ( std::vector<data_t>
			    ( geometry::Geometry<DIM>(e).bnd_vol(i)*data_t::storage_size ),
			    std::vector<data_t>
			    ( geometry::Geometry<DIM>(e).bnd_vol(i)*data_t::storage_size ) ) );
	recv_buff_.push_back
	  (std::make_pair ( std::vector<data_t>
			    ( geometry::Geometry<DIM>(e).bnd_vol(i)*data_t::storage_size ),
			    std::vector<data_t>
			    ( geometry::Geometry<DIM>(e).bnd_vol(i)*data_t::storage_size ) ) );
      }
    }

    static void init(int argc, char *argv[]) {
      int rc = MPI_Init(&argc,&argv);
      if (rc != MPI_SUCCESS) {
	printf ("Error starting MPI program. Terminating.\n");
	MPI_Abort(MPI_COMM_WORLD, rc);
      }
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs_);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

      // initialize neighbors *** only z direction ***
      nb_[0].first  = rank_;
      nb_[0].second = rank_;
      nb_[1].first  = rank_;
      nb_[1].second = rank_;
      nb_[2].first  = rank_;
      nb_[2].second = rank_;
      nb_[3].first  = (rank_ == numprocs_ - 1) ? 0 : rank_ + 1;
      nb_[3].second = (!rank_) ? numprocs_ - 1 : rank_ - 1;
    };

    ~Communicator() {
      MPI_Finalize();
    }

    const int& rank() { return rank_; }
    const nb_t& nb() { return nb_; }
    const buffer_t& send_buff() { return send_buff_; }
    const buffer_t& recv_buff() { return recv_buff_; }

    // Performs communication
    void operator()() {
      MPI_Request* r;
      for( Direction mu(0); mu.is_good(); ++mu)
	if (nb[int(mu)].second != rank) {

	  MPI_Isend(&send_buff[int(mu)].first[0],
		    send_buff[int(mu)].first.size()*data_t::storage_size,
		    MPI_DOUBLE, nb_[int(mu)].first,
		    int(mu), MPI_COMM_WORLD, &r);

	  MPI_Isend(&send_buff[int(mu)].second[0],
		    send_buff[int(mu)].second.size()*data_t::storage_size,
		    MPI_DOUBLE, nb_[int(mu)].second,
		    int(mu)+DIM, MPI_COMM_WORLD, &r);
	  // std::cout << "rank " << rank
	  //	    << "sends data on direction " << int(mu)
	  //	    << "to rank " << nb[int(mu)].second
	  //	    << std::endl;
	}

      for( Direction mu(0); mu.is_good(); ++mu)
	if (nb[int(mu)].first != rank) {

	  MPI_Status status;
	  MPI_Recv(&recv_buff[int(mu)].first[0],
		   recv_buff[int(mu)].first.size()*data_t::storage_size,
		   MPI_DOUBLE, nb_[int(mu)].first, int(mu)+DIM, MPI_COMM_WORLD,
		   &status);

	  MPI_Recv(&recv_buff[int(mu)].second[0],
		   recv_buff[int(mu)].second.size()*data_t::storage_size,
		   MPI_DOUBLE, nb_[int(mu)].second, int(mu), MPI_COMM_WORLD,
		   &status);

	  // std::cout << "rank " << rank
	  //	    << "receives data along direction " << int(mu)
	  //	    << "from rank " << nb[int(mu)].second
	  //	    << std::endl;
	}
    }

    // Check if any communication is still on its way
    bool test() { return false; }

  };

}

#endif //COMMUNICATOR_H
