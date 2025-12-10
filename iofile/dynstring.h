
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
#ifndef STRING_MK_H
#define STRING_MK_H

#ifdef _MSC_VER
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#endif

#include<fstream>
#include<iostream>
#include<string>
#include <cstring> //pour strlen



#define DIM 255
using namespace std;

class dynstring {
public:
  //Appel du constructeur avec la taille par défaut
    dynstring(const unsigned int Size=DIM){
        sz = Size;
        len = 0;
        list = new char[sz+1];
    }

  //constructeur avec un nombre de caractères de longeur L
  //la taille du string et donnée par Size;
    dynstring(unsigned int Length, char car, unsigned int Size=0) {
        if (!Size) Size = Length;
        sz = Size;
        if (Length >= sz) Length = sz;
        len = Length;
        list = new char [sz+1];
        for (char * Pch=list ; Length--; *Pch++ = car)
            *Pch = 0;
    }
  //constructeur avec une chaine de caractères
  dynstring(char * Ch, unsigned int Size=0)
    {
      if (!Size) Size = strlen(Ch);
      sz = Size;
      list = new char [sz+1];
      *this = Ch;
    }

  dynstring(const char * Ch)
  {
	  unsigned int Size = strlen(Ch);
	  sz = Size;
	  list = new char[sz + 1];
	  *this = Ch;
  }

    //destructeur
    virtual ~dynstring(){delete [] list;}


    //mÈthodes diverses
    //----------------------------------------
    int length() const {return len;}		//retourne la longeur de la chaine
    int size()  const {return sz;}		    //retourne la taille max de la chaine
    void empty() {len = 0; *list = 0;}	    //retourne la taille max de la chaine
    char * liste()  const {return list;}	//retourne la chaine de caractère

    // retourne 1 si les deux portions de
    //string sont identiques, sinon 0
    //----------------------------------------
    int equal(  const dynstring & b,
                unsigned int dest, unsigned int src, unsigned int Long) const {
        for(unsigned int i=0; i<Long; i++)
        if(list[dest+i] != b.list[src+i]) return 0;
        return 1;
    }

  //----------------------------------------
    dynstring * clone()			//retourne un ptr sur un clone du string
    {
        dynstring * Pstr = new dynstring (*this);
        Pstr->list = new char [sz+1];
        for (unsigned int i=0; i<sz; i++)
           Pstr->list[i] = list[i];
        return Pstr;
        // return new dynstring(*this);
    }

    void copy (const dynstring & orig)  	//copie le string spécifié
    {
        unsigned int Size = orig.sz;
        if (!(sz > Size)) Size = sz;
        len = orig.len;
        if (!(len < Size)) len = Size;
        for(unsigned int i=0; i<Size; i++)
            list[i] = orig.list[i];
    }

    //operateurs +
    //----------------------------------------
    dynstring & operator + (char Car)
    {
     dynstring * Temp = new dynstring (1);
     *Temp = Car;
     return *this + Temp->list;
    }


  //----------------------------------------
    dynstring & operator + (char * Ch)
    {
      //dynstring Temp (list, DIM);
      dynstring * Temp = new dynstring (list, sz+strlen (Ch));
      strncat((char*)Temp->list, Ch, sz+1);
      Temp->len = strlen (Temp->list);
      return *Temp;
    }

    //----------------------------------------
    dynstring & operator + (const char * Ch)
    {
        //dynstring Temp (list, DIM);
        dynstring * Temp = new dynstring(list, sz + strlen(Ch));
        strncat((char*)Temp->list, Ch, sz + 1);
        Temp->len = strlen(Temp->list);
        return *Temp;
    }
    //----------------------------------------
    dynstring & operator + (dynstring & Str) {return *this + Str.list;}


    //operateurs =
    //----------------------------------------
    dynstring & operator = (char Car)
    {
        list[0] = Car;
        list[1] = 0;
        len = 1;
        return *this;
    }
    //----------------------------------------
    dynstring & operator = (const char * Ch)
    {
      strncpy (list, Ch, sz+1);
      len = strlen(list);
      return *this;
    }
    //----------------------------------------
    dynstring & operator = (const dynstring & Str)
    {
      copy (Str);
      return *this;
    }


    //operateurs []
    char & operator [] (unsigned int Idx)
    {
      if (Idx > len) Idx = len+1;
      return list[Idx];
    }

    //Surcharge des operateurs d'extraction et d'inclusion
    //----------------------------------------
    friend ostream & operator << (ostream & Out, dynstring & Str) {
        return Out << Str.list;
    }
    //----------------------------------------
    friend istream & operator >> (istream & In, dynstring & Str){
        In >> Str.list; Str.len = strlen (Str.list); return In;
    }

    // surcharge d'opérateur de comparaison
    friend int operator==(const dynstring & a, const dynstring & b){
        //if((a.size() != b.size())||(a.length() != b.length())) return 0;
        if(a.length() != b.length()) return 0;
        return a.equal(b,0,0,a.length());
    }

    // surcharge d'opérateur de différence
    friend int operator!=(const dynstring & a, const dynstring & b){
        //if((a.length() != b.length()) || (a.size() != b.size())) return 1;
        if (a.length() != b.length()) return 1;
        return a.equal(b,0,0,a.length()) ? 0 : 1;
    }

    //------------------------------------------
protected:
  unsigned int len;		// number of elements
  unsigned int sz;		// dimensions
  char * list;			// character list
};


#endif
