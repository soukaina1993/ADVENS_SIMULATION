
/* ----------------------------------------------------------------------------
 DESCRIPTION:
  This defines the array base class.  Be careful what you put into one of these
arrays!  This class can be used only on objects that have:

  a default constructor (takes no arguments)
  operator=
  operator==
  operator!=

In other words, use only primitive objects in this array (e.g. int, float,
pointers, etc)

The constructor will try to initialize to zero, but only if the type is right.

We don't do any over-allocation, so resizing can be expensive.
No error checking on the copy, so don't walk over end of array!

---------------------------------------------------------------------------- */
#ifndef ARRAY_TPLT_H
#define ARRAY_TPLT_H

#ifdef _MSC_VER
    #pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#endif

template <class T>
class Array {
public:
  //Appel du constructeur avec la taille de l'objet
  Array(unsigned int s) : sz(s), a(sz ? new T[sz] : 0)
    {for(unsigned int i=0; i<sz; i++) a[i] = (T)0; }

  //Appel du constructeur avec un array sp�cifi�
  Array(const Array<T> & orig){sz=0; a=(T *)0; copy(orig);}

  //Surcharge de l'op�rateur d'�galit�
  Array<T> & operator=(const Array<T> & orig){copy(orig); return *this;}
  Array<T> & operator=(const T array []) // no err checks!
    {for(unsigned int i=0; i<sz; i++) a[i] = *(array+i); return *this;}

  //destructeur
  virtual ~Array(){delete [] a;}

  //retourne un ptr sur un clone de l'objet. L'utilisateur est
  //responsable de la m�moire allou�e pour l'appel de cette fonction
  Array<T> * clone(){return new Array<T>(*this);}
  operator const T * () const {return a;}
  operator T * () {return a;}

  //retourne une reference sur le i�me element de l'array
  const T & operator[](unsigned int i) const {return a[i];}
  T & operator[](unsigned int i) {return a[i];}
  const T & at(unsigned int i) const { return a[i]; }
  T & at(unsigned int i) { return a[i]; }

 //duplique l'array specifie
  void copy(const Array<T> & orig){
    size(orig.sz);
    for(unsigned int i=0; i<sz; i++)
      a[i] = orig.a[i];
  }

  //duplique une partie de l'array specifie. length elements a partir de
  //la position src sont copi� dans l'objet a partir de la position dest.
  //s'il y a pas assez d'espace dans l'objet, les �l�ments en exc�s
  //ne seront copies
  void copy(const Array<T> & orig, unsigned int dest,
	    unsigned int src, unsigned int length){
    for(unsigned int i=0; i<length; i++) a[dest+i] = orig.a[src+i];
  }

 //deplace le nombre d'elements specifies par length de la position src
 //a la position dest
  void move(unsigned int dest, unsigned int src, unsigned int length){
    if(src > dest)
      for(unsigned int i=0; i<length; i++) a[dest+i] = a[src+i];
    else if(src < dest)
      for(unsigned int i=length-1; i!=0; i--) a[dest+i] = a[src+i];
  }

 //�change les positions i et j dans l'objet
  void swap(unsigned int i, unsigned int j){T tmp=a[j]; a[j]=a[i]; a[i]=tmp;}

 //retourne la taille l'objet (l'array)
  unsigned int size() const {return sz;}

 //Sp�cifie la taille de l'objet (l'array)
  unsigned int size(unsigned int n){
    if(n == sz) return sz;
    T * tmp = (n ? new T[n] : 0);
    for(int i=((n < sz) ? n-1 : sz-1); i>=0; i--) tmp[i] = a[i];
    delete [] a;
    a = tmp;
    return sz=n;
  }

  // retourne 1 si les deux portions d'arrays sont identiques, sinon 0
  int equal(const Array<T> & b,
	    unsigned int dest, unsigned int src, unsigned int length) const {
    for(unsigned int i=0; i<length; i++)
      if(a[dest+i] != b.a[src+i]) return 0;
    return 1;
  }

  //Mettre la valeur val � la position pos et retourne un 1.
  //Retourne un 0 si la position sp�cifi�e > � la taille de l'objet
  int Set(unsigned int pos, T val )
  {
	  if (pos < sz){
	  	a[pos]=val;
	  	return 1;
	  }else return 0;
  }

  //retourne la valeur du i�me �l�ments de l'array
  T & Get(unsigned int pos)
  {if (pos < sz) return a[pos]; else return a[0];}
 //------------------------------------------

protected:
  unsigned int sz;		// number of elements
  T * a;			// the contents of the array
};

//------------------------------------------
// surcharge d'op�rateur de comparaison
template <class T> int
operator==(const Array<T> & a, const Array<T> & b){
  if(a.size() != b.size()) return 0;
  return a.equal(b,0,0,a.sz);
}
//------------------------------------------
// surcharge d'op�rateur de diff�rence
template <class T> int
operator!=(const Array<T> & a, const Array<T> & b){
  if(a.size() != b.size()) return 1;
  return a.equal(b,0,0,a.sz) ? 0 : 1;
}

typedef Array<double>	ArrayOfDouble;
typedef Array<float>	ArrayOfFloat;
typedef Array<int>	ArrayOfInt;
typedef Array<char>	ArrayOfChar;

#endif
