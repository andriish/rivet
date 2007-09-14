// -*- C++ -*-
#ifndef RIVET_Cmp_HH
#define RIVET_Cmp_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Cmp.fhh"
#include <typeinfo>


namespace Rivet {

   /// Cmp is a helper class to be used when checking the ordering of two
   /// objects. When implicitly converted to an integer the value will be
   /// negative if the two objects used in the constructor are ordered and
   /// positive if they are not. Zero will be returned if they are equal.
   ///
   /// The main usage of the Cmp class is if several variables should be
   /// checked for ordering in which case several Cmp objects can be
   /// combined as follows: <code>cmp(a1, a2) || cmp(b1, b2) || cmp(c1,
   /// c2)</code> where cmp is a global function for easy creation of Cmp
   /// objects.
  template <typename T>
  class Cmp {
    
  public:
    
    /// Enumerate the possible states of a Cmp object.
    enum CmpState {
      undefined = -2,  //< Undefined state.
      ordered = -1,    //< The two corresponding objects are ordered.
      equivalent = 0,  //< The two corresponding objects are equivalent.
      unordered = 1    //< The two corresponding objects are unordered.
    };
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    inline Cmp(const T & t1, const T & t2);
    
    /// The copy constructor.
    template <typename U>
    inline Cmp(const Cmp<U> &);
    
    /// The destructor is not virtual since this is not intended to be a base class.
    inline ~Cmp() {};
    
    /// The assignment operator.
    template <typename U>
    inline const Cmp<T> & operator=(const Cmp<U> &);
    
    //@}
    
  public:
    
    /// Automatically convert to an integer. 
    inline operator int () const;
    
    /// If this state is equivalent, set this state to the state of \a c.
    template <typename U>
    inline const Cmp<T> & operator||(const Cmp<U> & c) const;
    
  private:
    
    /// Perform the actual comparison if necessary.
    inline void compare() const;
    
  private:
    
    /// The state of this object.
    mutable int value;
    
    /// The objects to be compared.
    pair<const T *, const T *> objects;
    
  };
  
  /// Global helper function for easy creation of Cmp objects.
  template <typename T>
  inline Cmp<T> cmp(const T & t1, const T & t2) {
    return Cmp<T>(t1, t2);
  }
  
   /// Specialization of the Cmp helper class to be used when checking the
   /// ordering of two Projection objects. When implicitly converted to an
   /// integer the value will be negative if the two objects used in the
   /// constructor are ordered and positive if they are not. Zero will be
   /// returned if they are equal. This specialization uses directly the
   /// virtual compare() function in the Projection class.
   ///
   /// The main usage of the Cmp class is if several variables should be
   /// checked for ordering in which case several Cmp objects can be
   /// combined as follows: <code>cmp(a1, a2) || cmp(b1, b2) || cmp(c1,
   /// c2)</code> where cmp is a global function for easy creation of Cmp
   /// objects.
  template <>
  class Cmp<Projection> {
    
  public:
    
    /// Enumerate the possible states of a Cmp object.
    enum CmpState {
      undefined = -2, //< Undefined state.
      ordered = -1,   //< The two corresponding objects are ordered.
      equivalent = 0,  //< The two corresponding objects are equivalent.
      unordered = 1   //< The two corresponding objects are unordered.
    };
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    inline Cmp(const Projection & p1, const Projection & p2);
    
    /// The copy constructor.
    template <typename U>
    inline Cmp(const Cmp<U> &);
    
     /// The destructor is not virtual since this is not intended to be a base class.
    inline ~Cmp() {};
    
     /// The assignment operator.
    template <typename U>
    inline const Cmp<Projection> & operator=(const Cmp<U> &);
    
    //@}
    
  public:
    
    /// Automatically convert to an integer. 
    inline operator int () const;
    
    /// If this state is equaivalent, set this state to the state of \a c.
    template <typename U>
    inline const Cmp<Projection> & operator||(const Cmp<U> & c) const;
    
  private:
    
    /// Perform the actual comparison if necessary.
    inline void compare() const;
    
  private:
    
    /// The state of this object.
    mutable int value;
    
    /// The objects to be compared.
    pair<const Projection *, const Projection *> objects;
    
  };
  


  /////////////////////////////////////////////////////////



   /// Global helper function for easy creation of Cmp objects.
  inline Cmp<Projection> pcmp(const Projection & p1, const Projection & p2) {
    return Cmp<Projection>(p1, p2);
  }
  
  
  template <typename T>
  inline Cmp<T>::Cmp(const T & t1, const T & t2)
    : value(undefined), objects(&t1, &t2) {}
  
  template <typename T>
  template <typename U>
  inline Cmp<T>::Cmp(const Cmp<U> & x)
    : value(x.operator int()), objects(0, 0) {}
  
  template <typename T>
  template <typename U>
  inline const Cmp<T> & Cmp<T>::operator=(const Cmp<U> & x) {
    value = x.operator int();
    return *this;
  }
  
  template <typename T>
  template <typename U>
  inline const Cmp<T> & Cmp<T>::operator||(const Cmp<U> & x) const {
    compare();
    if ( value == equivalent ) value = x.operator int();
    return *this;
  }
  
  template <typename T>
  inline Cmp<T>::operator int() const {
    compare();
    return value;
  }
  
  template <typename T>
  inline void Cmp<T>::compare() const {
    if ( value == undefined ) {
      less<T> l;
      if ( l(*objects.first, *objects.second) ) value = ordered;
      else if ( l(*objects.second, *objects.first) ) value = unordered;
      else value = equivalent;
    }
  }
  
  inline Cmp<Projection>::Cmp(const Projection & p1, const Projection & p2)
    : value(undefined), objects(&p1, &p2) {}
  
  template <typename U>
  inline Cmp<Projection>::Cmp(const Cmp<U> & x)
    : value(x.operator int()), objects(0, 0) {}
  
  template <typename U>
  inline const Cmp<Projection> & Cmp<Projection>::operator=(const Cmp<U> & x) {
    value = x.operator int();
    return *this;
  }
  
  template <typename U>
  inline const Cmp<Projection> &
  Cmp<Projection>::operator||(const Cmp<U> & x) const {
    compare();
    if ( value == equivalent ) value = x.operator int();
    return *this;
  }
  
  inline Cmp<Projection>::operator int() const {
    compare();
    return value;
  }
  
  inline void Cmp<Projection>::compare() const {
    if ( value == undefined ) {
      const std::type_info & id1 = typeid(*objects.first);
      const std::type_info & id2 = typeid(*objects.second);
      if ( id1.before(id2) ) value = ordered;
      else if ( id2.before(id1) ) value = unordered;
      else {
        int c = objects.first->compare(*objects.second);
        if ( c < 0 ) value = ordered;
        else if ( c > 0 ) value = unordered;
        else value = equivalent;
      }
    }
  }
  
   
}


#endif
