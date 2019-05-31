#ifndef EGS_ARRAY2D_H_SEEN
#define EGS_ARRAY2D_H_SEEN

#include <stddef.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>

using namespace std;

/// \class  Array2D
/// \brief  Stores data of any type T in a 2D array
///
/// This class also provides a Fortran-like access operator ()
template <typename T>
class Array2D {
  public:
    /// \brief Default constructor, which does not allocate any memory
    Array2D(): xsize_(0), ysize_(0) {};

    /// \brief Constructor that allocates the memory
    Array2D(const size_t& nx, const size_t& ny):  xsize_(nx), ysize_(ny) {
        data_.resize(xsize_*ysize_);
    }

    /// \brief Constructor that allocates and initializes the data to a constant t
    Array2D(const size_t& nx, const size_t& ny, const T& t):  xsize_(nx), ysize_(ny) {
        data_.resize(xsize_*ysize_ , t);
    }

      /// \brief Destructor that frees up the memory
    ~Array2D() {data_.clear();}

    /// \brief Function to clear the memory
    void Clear() {
      xsize_ = 0;
      ysize_ = 0;
      data_.clear();
    }

    /// \brief Returns size in the x-direction
    size_t XSize() const {return xsize_;}
    /// \brief Returns size in the y-direction
    size_t YSize() const {return ysize_;}

    /// \brief Resizes the array
    /// \warning In its current implementation, most of the original data
    /// will get lost if the xsize changes as this changes the indexing for all entries.
    void Resize(const size_t& nx, const size_t& ny) {
      xsize_ = nx;
      ysize_ = ny;
      data_.resize(xsize_*ysize_);
    }

    /// \brief Resizes the array and sets ALL entries to the specified value
    void Resize(const size_t& nx, const size_t& ny, const T& t) {
      data_.clear();
      xsize_ = nx;
      ysize_ = ny;
      data_.resize(xsize_*ysize_, t);
    }

    /// \brief Set all values in the array to the given value
    void SetValue(const T& t){
      for(size_t i=0; i < data_.size(); i++){
        data_[i] = t;
      }
    }

    /// \brief Return a pointer to the first element of the data in the
    /// vector so we can use it access the data in array format
    T* GetArrayPointer() {
      return &(data_[0]);
    }

    /// \brief Return a cont point to the first element of the data in the
    /// vector so we can use it access the data in array format
    const T* GetConstArrayPointer() const {
      return &(data_[0]);
    }

    /// \brief Fortran-like () operator to access values in the 2D data array
    T& operator()(size_t ix,size_t iy) {return data_[ix+xsize_*iy];}

    /// \brief Fortran-like () const operator to access values in the 2D data array
    const T& operator()(size_t ix,size_t iy) const {return data_[ix+xsize_*iy];}

  private:

    /// \brief Copy constructor, which is made private so it would not be used inadvertently
    /// (until we define a proper copy constructor)
    Array2D(const Array2D &obj) {};

    /// \brief Number of elements in the x-dimension
    size_t xsize_;
    /// \brief Number of elements in the y-dimension
    size_t ysize_;

    vector<T> data_;
};

#endif /* ARRAY2D_H_SEEN */
