/* This file contains a program to take all the datasets in each of
 * several HDF5 files and stick them all into one HDF5 file. This
 * might be useful if you have run a FM several times in limited
 * configurations and want to use teh resulting radiances for one
 * humongous retrieval. Fancy things like links and user-defined types
 * are not handled yet.

 * Compile it like this:
   cc -g -o h5cat h5cat.c -lhdf5 
 * and run it like this:
   ./h5cat [-v] [-a attno] infile1 [infile2] [infile3.....]  outfile
   
 * -v turns on the verbosity -- without this the program is quite silent
 * -a nattrs limits the number of attributes for each dataset that are
 *    passed on to nattrs. (This option only exists to get round an 
 *    occasional infinite loop bug.)

  *  Copyright (c) © H. C. Pumphrey 2002 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <unistd.h> /* required for use of getopt() */

/* Global as I can't be arsed to pass them to the iterator */
/* hid_t outfile_id; */
char **global_varlist;
char **global_grouplist;
int n_vars,n_groups;
int global_nattrs, global_attrcnt,verbose;

/* Function prototype for operator (function is later on in this file) */
herr_t operator(hid_t group_id, const char *member_name, void
		       *operator_data);
/* Function prototype for attribute operator (function is later on in
   this file) */
herr_t attr_op(hid_t loc_id, const char *attribute_name, void
		       *operator_data);

typedef struct opdat_t{
  int getall;
  hid_t newloc_id;
}opdat_t;

/*-----Main---------Main---------Main---------Main---------Main-------*/
int main(int argc , char *argv[] ){

  /*--------Local vars--------------------------*/
  char infile[80], outfile[80];
  hid_t file_id,outfile_id,dset_id,dspace_id;
  int error, nmembers,*idx,obj_type, stringlen, i,j,filei,nfiles;
  char obj_name[80],rootname[80];
  char *varstring, *dummy,*groupstring;
  char *invarstring,*ingroupstring;
  
  void *operator_data;
  opdat_t opdat; 
  hid_t *newloc_idp,newloc_id;
  /* Stuff for getopt */
  char optstring[]="vd:g:a:";
  int optval;
  char optvalc;
  extern char *optarg;
  extern int optind,opterr,optopt;

  /* end of getopt vars */
 
  /*------------------Executable----------------------------*/
  strcpy(rootname,"/\0");
  if (argc < 3){
    printf("Usage: h5cat [-v] [-a nattr] infile1 [infile2...]  outfile \n");
    printf("nattr is the number of attributes to copy for each dataset\n");
    printf(" -v : verbose \n");
    exit(3);
  }
  
  /* Loop through options */
  invarstring=NULL;
  ingroupstring=NULL;
  optval=0;
  verbose=0;
  nfiles=0;
  while(optval >= 0){
    optval=getopt(argc,argv,optstring);
    optvalc=(char)optval;
    if(optvalc=='v'){
      verbose=1;
      printf("Ohhh! Verbosity on! Now I can spout lots of rubbish at you\n");
    }else if(optvalc=='a'){
      sscanf(optarg,"%d",&global_nattrs);
      if(verbose)
	printf("Doing %d attrs for each dataset \n",global_nattrs);
    }
  }/* End of getopt loop */

  nfiles=argc - optind;
  strcpy(outfile,argv[argc-1]);
  
  if(verbose){
    printf("copying from ");
    for(j=optind;j<argc-1;j++){
      printf(" %s ",argv[j]);
    }
    printf(" to %s \n ",outfile);
  }

  if (verbose)printf("opening outfile %s: id  was", outfile);
  outfile_id=H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  if (verbose)printf("%d\n",outfile_id);

  /* Loop over input files */
  for (filei=0;filei < nfiles-1; filei++){
    strcpy(infile,argv[optind+filei]);
    if(verbose) printf("opening file %d: %s id  was..",filei,infile);
    file_id=H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT );
    if(verbose)printf("...%d\n",file_id);

    /* iterate over all things in root group */
    idx=NULL;
    opdat.newloc_id=outfile_id;
    opdat.getall=0;
    error=H5Giterate(file_id, rootname, idx, 
		     operator, (void *) &opdat ) ;
  
    error=H5Fclose(file_id ) ;
  }
  error=H5Fclose(outfile_id ) ;
  
} /* End of main() */

/*---End-of-Main------End-of-Main------End-of-Main------End-of-Main---*/


/* This is the thing that gets done to each item in the root group. */
herr_t operator(hid_t group_id, const char *member_name, void
		*operator_data){
  /* Local variables */
  H5G_stat_t statbuf, *statbufp; 
  hbool_t follow_link;
  hid_t *newloc_idp, newloc_id;
  int error, is_wanted,i;
  opdat_t opdat, *opdatp, newopdat;
  /* executable */
  opdatp= (opdat_t *) operator_data;
 
  newloc_id=opdatp->newloc_id;

  if(verbose) printf("operating on %s ... ",member_name);
    
  follow_link=0;
  statbufp=&statbuf;

  /* Call to HDF5 to find what type of item this is. The results are returned
   * in the structure statbuf, pointed to by statbufp */
  error=H5Gget_objinfo(group_id, member_name, follow_link, statbufp );
  
  /* Main if block  to process each of the possible item types */
  if(statbuf.type ==  H5G_GROUP){
    hid_t oldsubgroup_id, newsubgroup_id, *newsubgroup_idp;
    int *idx;
    opdat_t newsubgroup;
    size_t size_hint;
    if(verbose) printf("It's a group: Looking into it\n");
    /*oldsubgroup_id=H5Gopen(group_id,member_name );*/
    size_hint=40;
    newsubgroup.newloc_id=H5Gcreate(newloc_id, member_name, size_hint );
    newsubgroup.getall=newopdat.getall;
    idx=NULL;
    error=H5Giterate(group_id, member_name, idx, 
		     operator, (void *) &newsubgroup ) ;

    /*error=H5Gclose(oldsubgroup_id); */
    error=H5Gclose(newsubgroup.newloc_id); 
    if(verbose) printf("Finished with group %s\n",member_name);
  }
  else if(statbuf.type ==  H5G_DATASET){
    /* OK, it's a dataset. Now we need to process it to write it out */
    /* variables local to this block (?is this standard C? It seems to work.)*/
    hid_t old_dspace, new_dspace, old_dset, new_dset, type_id ,plist_id;
    void *databuf;
    int n_elements;
    int *aidx;
    void *aopdata;
    hid_t *new_dsetp;
    /* executable for this bock */
    if(verbose) printf("It's a dataset: \n");
    old_dset=H5Dopen(group_id, member_name );
    /*printf("Opened dataset with id=%d\n",old_dset); */
    old_dspace=H5Dget_space(old_dset) ;
    /*printf("Opened dataspace with id=%d\n",old_dspace);*/
    plist_id=H5Dget_create_plist(old_dset); 
    new_dspace=H5Scopy(old_dspace);
    type_id=H5Dget_type(old_dset); 
    
    /* This next line makes a new dataset in the new file with an
     * identical dataspace to the old dataset in the old file */
    /* needs fixing to work for extendible datasets */
    new_dset=H5Dcreate(newloc_id, member_name, type_id, 
		       new_dspace, plist_id );

    /* Now we have to move the actual data. To do that, we have to malloc
       a suitably humongous buffer */
    n_elements=H5Sget_simple_extent_npoints(new_dspace); 
    databuf=(void *) malloc(n_elements*H5Tget_size(type_id) );
    error=H5Dread(old_dset,type_id, H5S_ALL,H5S_ALL ,
		  H5P_DEFAULT,databuf ); 
    error=H5Dwrite(new_dset,type_id,H5S_ALL,H5S_ALL, 
		   H5P_DEFAULT,databuf ) ;
    /* Process all attributes of this dataset */
    aidx=NULL;
    new_dsetp=&new_dset;
    aopdata=(void *) new_dsetp;
    global_attrcnt=0;
    if(global_nattrs > 0)
      error=H5Aiterate(old_dset, aidx, attr_op, aopdata ) ;


    free(databuf);
    error=H5Sclose(old_dspace );
    error=H5Sclose(new_dspace );
    error=H5Dclose(old_dset); 
    error=H5Dclose(new_dset); 
  }
  else if(statbuf.type ==  H5G_LINK){
    if(verbose) printf("It's a link: ignore\n");
  }
  else if(statbuf.type ==  H5G_TYPE){
    if(verbose) printf("It's a type: ignore\n");
  }
  else {
    printf("I really don't know what the fsck %s is: it has type:%d\n",
	   member_name,statbuf.type);
  }
  /* End of main if block  to process each of the possible item types */  

  return 0;

}

herr_t attr_op(hid_t loc_id, const char *attr_name, void *operator_data){

  hid_t attr_id, new_attr_id,new_dset_id, *new_dset_idp;
  hid_t type_id, dspace_id,new_space,space_id,new_dspace;
  int error,n_elements;
  void *databuf;
  new_dset_idp=(hid_t *) operator_data;
  new_dset_id= *new_dset_idp;


  printf("Processing attribute %s\n",attr_name);
  
  attr_id=H5Aopen_name(loc_id, attr_name ) ;
  type_id=H5Aget_type(attr_id);
  space_id=H5Aget_space(attr_id);
  new_dspace=H5Scopy(space_id);
  new_attr_id=H5Acreate(new_dset_id,attr_name,type_id, space_id, H5P_DEFAULT);
  
  n_elements=H5Sget_simple_extent_npoints(new_dspace); 
  databuf=(void *) malloc(n_elements*H5Tget_size(type_id) );
  error=H5Aread(attr_id,type_id, databuf ); 
  error=H5Awrite(new_attr_id,type_id, databuf ); 

  
  error=H5Aclose(attr_id);


  global_attrcnt++;
  if(global_attrcnt < global_nattrs){
    return 0; /* next attribute will be processed */
  }else{
    return 1; /* No more attributes will be processed */
  }
}
