/* 
 * Open-boundary injection code which appears to derive from example deck in
 * ./sample/reconnection/open-collisional/head/open-collisional
 * Reference: Daughton, Scudder, Karimabadi (2006, PoP).
 */

//  Safe allocate
#define ALLOCATE(A,LEN,TYPE)                                            \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));

#define DUMP_INJECTOR(side,dim1,dim2,set)                               \
  {                                                                     \
    FileIO fileIO;                                                      \
    FileIOStatus status;                                                \
    char filename[40];                                                  \
    sprintf(filename, "restart" #set "/" #side ".%d", int(rank()));     \
    status= fileIO.open(filename, io_write);                            \
    size_t len = dim1*dim2*nsp;                                         \
    fileIO.write(global->c##side,len);                                  \
    fileIO.write(global->n##side,len);                                  \
    fileIO.write(global->u##side,3*len);                                \
    fileIO.write(global->p##side,9*len);                                \
    fileIO.close();                                                     \
  }

#define READ_INJECTOR(side,dim1,dim2,set)                               \
  {                                                                     \
    FileIO fileIO;                                                      \
    FileIOStatus status;                                                \
    char filename[40];                                                  \
    sprintf(filename, "restart" #set "/" #side ".%d",int(rank()));      \
    status= fileIO.open(filename, io_read);                             \
    size_t len = dim1*dim2*nsp;                                         \
    fileIO.read(global->c##side,len);                                   \
    fileIO.read(global->n##side,len);                                   \
    fileIO.read(global->u##side,3*len);                                 \
    fileIO.read(global->p##side,9*len);                                 \
    fileIO.close();                                                     \
  }

#define DEFINE_INJECTOR(side,dim1,dim2)         \
  {                                             \
    size_t len = dim1*dim2*nsp;                 \
    ALLOCATE(global->c##side,len,int);          \
    ALLOCATE(global->n##side,len,double);       \
    ALLOCATE(global->u##side,3*len,double)      \
    ALLOCATE(global->p##side,9*len,double);			\
  }

# define DUMP_INJECTORS(set)                            \
  if (global->center) DUMP_INJECTOR(center,ny,nz,set);
//  if (global->left) DUMP_INJECTOR(left,ny,nz,set);      \
//  if (global->right) DUMP_INJECTOR(right,ny,nz,set);    \
//  if (global->top) DUMP_INJECTOR(top,ny,nx,set);        \
//  if (global->bottom) DUMP_INJECTOR(bot,ny,nx,set);
      
// Define Fortran style indexing of arrays
#define INDEX_FORTRAN_4(a,b,c,d,al,ah,bl,bh,cl,ch,dl,dh)                \
  ((a)-(al) + ((ah)-(al)+1)*((b)-(bl) +  ((bh)-(bl)+1)*(((c)-(cl)) + ((ch)-(cl)+1)*((d)-(dl)))))
#define INDEX_FORTRAN_5(a,b,c,d,e,al,ah,bl,bh,cl,ch,dl,dh,el,eh)        \
  ((a)-(al) + ((ah)-(al)+1)*((b)-(bl) +  ((bh)-(bl)+1)*(((c)-(cl)) + ((ch)-(cl)+1)*(((d)-(dl)) + ((dh)-(dl)+1)*((e)-(el))))))

#define ccenter(n,i,j) global->ccenter[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define ncenter(n,i,j) global->ncenter[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]   
#define ucenter(a,n,i,j) global->ucenter[INDEX_FORTRAN_4(a,n,i,j,1,3,1,nsp,1,nz,1,ny)]
#define pcenter(a,b,n,i,j) global->pcenter[INDEX_FORTRAN_5(a,b,n,i,j,1,3,1,3,1,nsp,1,nz,1,ny)]

#define cleft(n,i,j) global->cleft[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define nleft(n,i,j) global->nleft[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]   
#define uleft(a,n,i,j) global->uleft[INDEX_FORTRAN_4(a,n,i,j,1,3,1,nsp,1,nz,1,ny)]
#define pleft(a,b,n,i,j) global->pleft[INDEX_FORTRAN_5(a,b,n,i,j,1,3,1,3,1,nsp,1,nz,1,ny)]

#define cright(n,i,j) global->cright[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]
#define nright(n,i,j) global->nright[INDEX_FORTRAN_3(n,i,j,1,nsp,1,nz,1,ny)]   
#define uright(a,n,i,j) global->uright[INDEX_FORTRAN_4(a,n,i,j,1,3,1,nsp,1,nz,1,ny)]
#define pright(a,b,n,i,j) global->pright[INDEX_FORTRAN_5(a,b,n,i,j,1,3,1,3,1,nsp,1,nz,1,ny)]

#define vth(n) global->vth[INDEX_FORTRAN_1(n,1,nsp)]
#define q(n) global->q[INDEX_FORTRAN_1(n,1,nsp)]
#define npcenter(n) global->npcenter[INDEX_FORTRAN_1(n,1,nsp)]
#define npleft(n) global->npleft[INDEX_FORTRAN_1(n,1,nsp)]
#define npright(n) global->npright[INDEX_FORTRAN_1(n,1,nsp)]
