
/*#*******************	
*  Adjacency Members *
**********************/


/** Return a boolean that indicate if the face is complex.
    @param j Index of the edge
	@return true se la faccia e' manifold, false altrimenti
*/
inline bool IsManifold( const int j ) const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	return ( F(j)==this || this == F(j)->F(Z(j)) );
#endif
	return true;
	assert(0);
}

/** Return a boolean that indicate if the j-th edge of the face is a border.
	@param j Index of the edge
	@return true if j is an edge of border, false otherwise
*/
inline bool IsBorder( const int j ) const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	return F(j)==this;
#endif
	return true;
	assert(0);
}


/// This function counts the boreders of the face
inline int BorderCount() const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	int t = 0;
	if( IsBorder(0) ) ++t;
	if( IsBorder(1) ) ++t;
	if( IsBorder(2) ) ++t;
	return t;
#endif
	assert(0);
	return 3;
}


/// This function counts the number of incident faces in a complex edge
inline int ComplexSize(const int e) const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	int cnt=0;
	FACE_TYPE *fi=(FACE_TYPE *)this;
	int nzi,zi=e;
	do
	{
		nzi=fi->Z(zi);
		fi=fi->F(zi);
		zi=nzi;
		++cnt;
	}
	while(fi!=this);
	return cnt;
#endif
	assert(0);
	return 2;
}

/*Funzione di detach che scollega una faccia da un ciclo 
(eventualmente costituito da due soli elementi) incidente su un edge*/
/** This function detach the face from the adjacent face via the edge e. It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param e Index of the edge
*/
void Detach(const int e)
{
	typedef FEdgePosB< FACE_TYPE > ETYPE;

	assert(!IsBorder(e));
	ETYPE EPB(this,e);  // la faccia dall'altra parte
	EPB.NextF();
	int cnt=0;
	while ( EPB.f->F(EPB.z) != this)
	{ 
		assert(!IsManifold(e));   // Si entra in questo loop solo se siamo in una situazione non manifold.
		assert(!EPB.f->IsBorder(EPB.z));
		EPB.NextF();
		cnt++;
	}
	assert(EPB.f->F(EPB.z)==this);

	EPB.f->F(EPB.z) = F(e);
	EPB.f->Z(EPB.z) = Z(e);
	
	F(e) = this;
	Z(e) = e;

	EPB.f->SetM();
	this->SetM();
}

void OldDetach(const int e)
{
	typedef EdgePosB< FACE_TYPE > ETYPE;

	assert(!IsBorder(e));
	ETYPE EPB(this,e);
	ETYPE TEPB(0,-1);
	EPB.NextF();
	while ( EPB.f != this)
	{
		TEPB = EPB;
		assert(!EPB.f->IsBorder(EPB.z));
		EPB.NextF();
	}
	assert(TEPB.f->F(TEPB.z)==this);
	TEPB.f->F(TEPB.z) = F(e);
	TEPB.f->Z(TEPB.z) = Z(e);
	F(e) = this;
	Z(e) = e;
	TEPB.f->SetM();
	this->SetM();
}


/** This function attach the face (via the edge z1) to another face (via the edge z2). It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param z1 Index of the edge
		@param f2 Pointer to the face
		@param z2 The edge of the face f2 
*/
void Attach(int z1, face_base *&f2, int z2)
{
	typedef FEdgePosB< FACE_TYPE > ETYPE;
	ETYPE EPB(f2,z2);
	ETYPE TEPB;
	TEPB = EPB;
	EPB.NextF();
	while( EPB.f != f2)  //Alla fine del ciclo TEPB contiene la faccia che precede f2
	{
		TEPB = EPB;
		EPB.NextF();
	}
	//Salvo i dati di f1 prima di sovrascrivere
	face_base *f1prec = this->F(z1);  
	int z1prec = this->Z(z1);
	//Aggiorno f1
	this->F(z1) = TEPB.f->F(TEPB.z);  
	this->Z(z1) = TEPB.f->Z(TEPB.z);
	//Aggiorno la faccia che precede f2
	TEPB.f->F(TEPB.z) = f1prec;
	TEPB.f->Z(TEPB.z) = z1prec;
}


void AssertAdj()
{
	assert(F(0)->F(Z(0))==this);
	assert(F(1)->F(Z(1))==this);
	assert(F(2)->F(Z(2))==this);

	assert(F(0)->Z(Z(0))==0);
	assert(F(1)->Z(Z(1))==1);
	assert(F(2)->Z(Z(2))==2); 
}
// Funzione di supporto
inline void Nexts( face_base *&f,int &z )
{
    int t;
    t = z;
    z = (*f).Z(z);
    f = (*f).F(t);
}

/** This function change the orientation of the face. Inverting the index of two vertex 
@param z Index of the edge
*/
void Swap ( const int z )
{

  int i;
  face_base *tmp, *prec;
  int t, precz;

  swap ( V((z  )%3),V((z+1)%3));

  if( OBJ_TYPE & (OBJ_TYPE_A|OBJ_TYPE_S ) )
  {
	swap ( F((z+1)%3),F((z+2)%3));
	swap ( Z((z+1)%3),Z((z+2)%3));

	for(i = 1; i < 3; i++)
	{

      tmp = this;
      t = (z+i)%3;
      do {
					prec = tmp;
					precz = t;
					Nexts(tmp,t);
      }
      while (tmp != this);
  
      (*prec).Z(precz) = (z+i)%3;
    }
  }
}



// Stacca la faccia corrente dalla catena di facce incidenti sul vertice z, 
// NOTA funziona SOLO per la topologia VF!!!
// usata nelle classi di collapse
void VFDetach(int z)
{
	
	if(V(z)->Fp()==this )
	{
		int fz = V(z)->Zp();
		V(z)->Fp() = (face_from_vert_type *) F(fz);
		V(z)->Zp() = Z(fz);
	}
	else
	{
			VEdgePosB<FACE_TYPE> x,y;

		x.f = V(z)->Fp();
		x.z = V(z)->Zp();

		for(;;)
		{
			y = x;
			x.NextF();
			assert(x.f!=0);
			if(x.f==this)
			{
				y.f->F(y.z) = F(z);
				y.f->Z(y.z) = Z(z);
				break;
			}
		}
	}
}





}	 // end namespace


#endif

