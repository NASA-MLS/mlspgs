/* This file contains a program to extract some of the datasets in an
 * HDF5 file and put them into another HDF5 file. This can be useful
 * if you have a file containing many datasets and you only want one
 * or two, particularly if you want to transport the file across a
 * slow link, fit it on a small hard disc, read it into a program that
 * barfs on large files and can't do the subsetting itself (e.g. R) etc. etc. 

 * Compile it like this:
   cc -g -o h5subset h5subset.c -lhdf5 
 * and run it like this:
   ./h5subset infile outfile varfoo,varbar,varbaz........

 * .....to extract variables varfoo,varbar,varbaz from infile and put them
 * in outfile 

 * varfoo,varbar etc can be groups as well as datasets. 
 
 *  Copyright (c) © H. C. Pumphrey 2001 */


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
  int error, nmembers,*idx,obj_type, stringlen, i;
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
    printf("Usage: h5subset -a nattr -d dsets -g groups <infile> <outfile>\n");
    printf("where dsets is a comma separated list of datasets to copy\n"); 
    printf("groups is a comma-sep. list of groups to copy all of \n");
    printf("nattr is the number of attributes to copy for each dataset\n");
    printf(" -v : verbose \n");
    exit(3);
  }
  
  /* Loop through options */
  invarstring=NULL;
  ingroupstring=NULL;
  optval=0;
  verbose=0;
  while(optval >= 0){
    optval=getopt(argc,argv,optstring);
    optvalc=(char)optval;
    if(optvalc=='v'){
      verbose=1;
      printf("Ohhh! Verbosity on! Now I can spout lots of rubbish at you\n");
    }
    else if(optvalc=='d'){
      stringlen=strlen(optarg);
      invarstring=(char *) malloc(stringlen+3);
      varstring=(char *) malloc(stringlen+3);
      strcpy(invarstring,optarg);
      strcpy(varstring,optarg);
    }else if(optvalc=='g'){
      if(verbose) printf("arg to -g was%s\n",optarg);
      stringlen=strlen(optarg);
      ingroupstring=(char *) malloc(stringlen+3);
      groupstring=(char *) malloc(stringlen+3);
      strcpy(groupstring,optarg); 
      strcpy(ingroupstring,optarg); 
    }else if(optvalc=='a'){
      sscanf(optarg,"%d",&global_nattrs);
      if(verbose)
	printf("Doing %d attrs for each dataset \n",global_nattrs);
    }
  }/* End of getopt loop */
  
  strcpy(infile,argv[argc-2]);
  strcpy(outfile,argv[argc-1]);

  if(verbose)printf("copying from %s to %s \n ",infile,outfile);
  /* split varstring into an array of variable names. */
  /* First, count the variables */
  n_vars=0;
  dummy=strtok(varstring,",");
  while(dummy !=NULL){
    n_vars++;
    dummy=strtok(NULL,",");
  }
  /* Now allocate an array of strings for them and put them in it */
  if(verbose)printf("That is %d variables\n ",n_vars);
  global_varlist=(char **) malloc(sizeof(char *)*n_vars);
  strcpy(varstring,invarstring);
  
  global_varlist[0]=strtok(varstring,",");
  if(verbose)printf("Will search for variables X%sX  ",global_varlist[0]);
  for(i=1;i<n_vars;i++){
    global_varlist[i]=strtok(NULL,",");
    if(verbose) printf(" X%sX ",global_varlist[i]);
  }
  if(verbose) printf("\n");

  /* Done making string array of  variable names */
  
  if(ingroupstring!=NULL){
    /* split groupstring into an array of variable names. */
    /* First, count the variables */
    n_groups=0;
    dummy=strtok(groupstring,",");
    while(dummy !=NULL){
      n_groups++;
      dummy=strtok(NULL,",");
    }
    /* Now allocate an array of strings for them and put them in it */
    
    if(verbose) printf("That is %d groups\n ",n_groups);
    global_grouplist=(char **) malloc(sizeof(char *)*n_groups);
    strcpy(groupstring,ingroupstring);
  
    global_grouplist[0]=strtok(groupstring,",");
    if(verbose) printf("Will search for groups X%sX  ",global_grouplist[0]);
    for(i=1;i<n_groups;i++){
      global_grouplist[i]=strtok(NULL,",");
      if(verbose) printf(" X%sX ",global_grouplist[i]);
    }
    if(verbose) printf("\n");
  }

  /* printf("opening file: %s id  was",infile);*/
  file_id=H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT );
  /*printf("%d\n",file_id);*/

  /*printf("opening outfile %s: id  was", outfile);*/
  outfile_id=H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  /*printf("%d\n",outfile_id);*/
  /* iterate over all things in root group */
  idx=NULL;
  opdat.newloc_id=outfile_id;
  opdat.getall=0;
  error=H5Giterate(file_id, rootname, idx, 
		   operator, (void *) &opdat ) ;
  
  error=H5Fclose(file_id ) ;
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
  if(opdatp->getall==1){
    if(verbose) printf("getting %s whatever\n",member_name);
    newopdat.getall=1;
    is_wanted=1;
  }
  else{
    if(verbose) printf("operating on %s ... ",member_name);
    is_wanted=0;
    newopdat.getall=0;
      for(i=0;i<n_vars;i++){
	if (strcmp(global_varlist[i],member_name)==0){
	  is_wanted=1;
	} 
      }
    for(i=0;i<n_groups;i++){
      if (strcmp(global_grouplist[i],member_name)==0){
	is_wanted=1;
	newopdat.getall=1;
      } 
    }
  }
  if(is_wanted==0){
    if(verbose) printf("not on wanted list\n");
  }
  else {
    
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
  }
  /* end of wanted-or-not if block */
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
