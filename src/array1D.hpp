#ifndef EGS_ARRAY1D_H_SEEN
#define EGS_ARRAY1D_H_SEEN

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>

using namespace std;

/// \class  Array1D
/// \brief  Stores data of any type T in a 1D array

template <typename T>
class Array1D {
  public:
    /// \brief Default constructor, which does not allocate any memory
    Array1D(): xsize_(0) {};

    /// \brief Constructor that allocates the memory
    Array1D(const size_t& nx): xsize_(nx) {
      data_.resize(xsize_);
    }

    /// \brief Constructor that allocates and initializes the data to a value t
    Array1D(const size_t& nx, const T& t): xsize_(nx) {
      data_.resize(xsize_, t);
    }

    /// \brief Assignment operator copies the data structure by value
    Array1D& operator=(const Array1D &obj) {
      xsize_ = obj.xsize_;
      data_ = obj.data_;
      return *this;
    }

    /// \brief Copy constructor
    Array1D(const Array1D &obj): xsize_(obj.xsize_), data_(obj.data_) {};

    /// \brief Destructor that frees up the memory
    ~Array1D() {data_.clear();}

    /// \brief Function to clear the memory
    void Clear() {
      xsize_ = 0;
      data_.clear();
    }

    /// \brief Returns size in the x-direction
    size_t XSize() const {return xsize_;}

    /// \brief Returns length 
    size_t Length() const {return xsize_;}

    /// \brief Resizes the array
    void Resize(const size_t& nx) {
      xsize_ = nx;
      data_.resize(xsize_);
    }

    /// \brief Resizes the array and sets ALL entries to the specified value
    /// \warning All original data will get lost if this function is used!
    void Resize(const size_t& nx, const T& t) {
      data_.clear();
      xsize_ = nx;
      data_.resize(xsize_, t);
    }

    /// \brief Set all values in the array to the given value
    void SetValue(const T& t){
      for(size_t i=0; i < data_.size(); i++){
        data_[i] = t;
      }
    }

    /// \brief Add element to the end of the vector
    void PushBack(const T& t){
       xsize_ += 1;
       data_.push_back(t);
    }

    /// \brief Return a pointer to the first element of the data in the
    /// vector so we can use it access the data in array format
    T* GetArrayPointer() {
      return &(data_[0]);
    }

    /// \brief Return a const point to the first element of the data in the
    /// vector so we can use it access the data in array format
    const T* GetConstArrayPointer() const {
      return &(data_[0]);
    }

    /// \brief Fortran-like () operator to access values in the 1D data array
    T& operator()(size_t ix) {return data_[ix];}

    /// \brief Fortran-like () const operator to access values in the 1D data array
    const T& operator()(size_t ix) const {return data_[ix];}

  private:

    /// \brief Number of elements
    size_t xsize_;

    /// \brief Data in the array with size = xsize_
    vector<T> data_;

};

#endif /* ARRAY1D_H_SEEN */
