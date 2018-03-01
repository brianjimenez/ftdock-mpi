/*
Copyright (C) 1997-2001 Imperial Cancer Research Technology
author: Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

*/

#include "structures.h"

int main( int argc , char *argv[] ) {

  /* index counters */

  int		i ;

  /* Command line options */

  char		*static_file_name ;
  char		*mobile_file_name ;
  char		*matrix_file_name ;
  char		*gf_type ;
  char		gf_one_letter_code ;

  char		*default_matrix_file_name ;

  /* File stuff */

  FILE		*matrix_file ;
  char		line_buffer[200] ;

  /* ftdock_*.dat file values */

  float		RPscore ;

  /* Structures */

  struct Structure	Static_Structure , Mobile_Structure ;

  /* Grid stuff */

  int		static_residue , static_atom ;
  int		mobile_residue , mobile_atom ;
  int		prox ;

  struct Matrix	Matrix ;

/************/

  /* Its nice to tell people what going on straight away */

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "\n          3D-Dock Suite (May 2000)\n" ) ;
  printf( "          Copyright (C) 1997-2001 Imperial Cancer Research Technology\n" ) ;
  printf( "          Biomolecular Modelling Laboratory\n" ) ;
  printf( "          Imperial Cancer Research Fund\n" ) ;
  printf( "          44 Lincoln's Inn Fields\n" ) ;
  printf( "          London WC2A 3PX\n" ) ;
  printf( "          +44 (0)20 7269 3565\n" ) ;
  printf( "          http://www.bmm.icnet.uk/\n\n" ) ;

  printf( "Starting RPDock residue level pair potential score program\n" ) ;

/************/

  /* Memory allocation */

  if( ( ( matrix_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( gf_type           = ( char * ) malloc (  10 * sizeof( char ) ) ) == NULL ) ||
      ( ( static_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( mobile_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ) {
    GENERAL_MEMORY_PROBLEM 
  }

/************/

  /* Command Line defaults */

  strcpy( matrix_file_name , "best.matrix" ) ;

  strcpy( gf_type , "G_DATA" ) ;
  gf_one_letter_code = 'g' ;

  default_matrix_file_name = "(default)" ;

  /* Command Line parse */

  for( i = 1 ; i < argc ; i ++ ) {

    if( strcmp( argv[i] , "-p1" ) == 0 ) {
      i ++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
        printf( "Bad command line\n" ) ;
        exit( EXIT_FAILURE ) ;
      }
      strcpy( static_file_name , argv[i] ) ;
    } else {
      if( strcmp( argv[i] , "-p2" ) == 0 ) {
        i ++ ;
        if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
          printf( "Bad command line\n" ) ;
          exit( EXIT_FAILURE ) ;
        }
        strcpy( mobile_file_name , argv[i] ) ;
      } else {
        if( strcmp( argv[i] , "-matrix" ) == 0 ) {
          i ++ ;
          if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
          }
          strcpy( matrix_file_name , argv[i] ) ;
          default_matrix_file_name = "(user defined)" ;
        } else {
          printf( "Bad command line\n" ) ;
          exit( EXIT_FAILURE ) ;
        }
      }
    }

  }

/************/

  /* Read in matrix */

  if( ( matrix_file = fopen( matrix_file_name , "r" ) ) == NULL ) {
    printf( "Could not open %s for reading.\nDying\n\n" , matrix_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }

  fgets( line_buffer , 199 , matrix_file ) ;
  if( strncmp( line_buffer , "Description" , 11 ) == 0 ) {
    sscanf( line_buffer , "Description       :: %s" , Matrix.description ) ;
  } else {
    printf( "Error in matrix file\nDying\n" ) ;
    exit( EXIT_FAILURE ) ;
  }
  fgets( line_buffer , 199 , matrix_file ) ;
  if( strncmp( line_buffer , "Distance cut off" , 16 ) == 0 ) {
    sscanf( line_buffer , "Distance cut off  :: %f" , &Matrix.distance ) ;
  } else {
    printf( "Error in matrix file\nDying\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

  for( i = 1 ; i <= 20 ; i ++ ) {

    fgets( line_buffer , 199 , matrix_file ) ;
    
    sscanf( line_buffer       , "%8f" , &Matrix.score[i][1] ) ;
    sscanf( line_buffer +   8 , "%8f" , &Matrix.score[i][2] ) ;
    sscanf( line_buffer +  16 , "%8f" , &Matrix.score[i][3] ) ;
    sscanf( line_buffer +  24 , "%8f" , &Matrix.score[i][4] ) ;
    sscanf( line_buffer +  32 , "%8f" , &Matrix.score[i][5] ) ;
    sscanf( line_buffer +  40 , "%8f" , &Matrix.score[i][6] ) ;
    sscanf( line_buffer +  48 , "%8f" , &Matrix.score[i][7] ) ;
    sscanf( line_buffer +  56 , "%8f" , &Matrix.score[i][8] ) ;
    sscanf( line_buffer +  64 , "%8f" , &Matrix.score[i][9] ) ;
    sscanf( line_buffer +  72 , "%8f" , &Matrix.score[i][10] ) ;
    sscanf( line_buffer +  80 , "%8f" , &Matrix.score[i][11] ) ;
    sscanf( line_buffer +  88 , "%8f" , &Matrix.score[i][12] ) ;
    sscanf( line_buffer +  96 , "%8f" , &Matrix.score[i][13] ) ;
    sscanf( line_buffer + 104 , "%8f" , &Matrix.score[i][14] ) ;
    sscanf( line_buffer + 112 , "%8f" , &Matrix.score[i][15] ) ;
    sscanf( line_buffer + 120 , "%8f" , &Matrix.score[i][16] ) ;
    sscanf( line_buffer + 128 , "%8f" , &Matrix.score[i][17] ) ;
    sscanf( line_buffer + 136 , "%8f" , &Matrix.score[i][18] ) ;
    sscanf( line_buffer + 144 , "%8f" , &Matrix.score[i][19] ) ;
    sscanf( line_buffer + 152 , "%8f" , &Matrix.score[i][20] ) ;

  }

  fclose( matrix_file ) ;
  
/************/

  /* Read in Structures from pdb files */
  Static_Structure = read_pdb_to_structure( static_file_name ) ;
  Mobile_Structure = read_pdb_to_structure( mobile_file_name ) ;

/************/

  /* Main program loop */

  /* No buffering at all */
  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  RPscore = 0 ;

  for( static_residue = 1 ; static_residue <= Static_Structure.length ; static_residue ++ ) {

    if( Static_Structure.Residue[static_residue].nc > 20 ) break ;

    for( mobile_residue = 1 ; mobile_residue <= Mobile_Structure.length ; mobile_residue ++ ) {

      if( Mobile_Structure.Residue[mobile_residue].nc > 20 ) break ;

      prox = 0 ;

      for( static_atom = 1 ; static_atom <= Static_Structure.Residue[static_residue].size ; static_atom ++ ) {
        for( mobile_atom = 1 ; mobile_atom <= Mobile_Structure.Residue[mobile_residue].size ; mobile_atom ++ ) {

          if( ( ( Static_Structure.Residue[static_residue].Atom[static_atom].coord[1] - Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[1] ) * ( Static_Structure.Residue[static_residue].Atom[static_atom].coord[1] - Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[1] ) ) + ( ( Static_Structure.Residue[static_residue].Atom[static_atom].coord[2] - Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[2] ) * ( Static_Structure.Residue[static_residue].Atom[static_atom].coord[2] - Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[2] ) ) + ( ( Static_Structure.Residue[static_residue].Atom[static_atom].coord[3] - Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[3] ) * ( Static_Structure.Residue[static_residue].Atom[static_atom].coord[3] - Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[3] ) ) < ( Matrix.distance * Matrix.distance ) ) prox = 1 ;

        }
      }

      if( prox == 1 ) RPscore += Matrix.score[Static_Structure.Residue[static_residue].nc][Mobile_Structure.Residue[mobile_residue].nc] ;

    }
  }

  fprintf( stdout, "%6s   %8.3f\n" , gf_type , RPscore ) ;

  printf( "\n" ) ;

/************/

  /* Free the memory */

  for( i = 1 ; i <= Static_Structure.length ; i ++ ) {
    free( Static_Structure.Residue[i].Atom ) ;
  }
  free( Static_Structure.Residue ) ;

  for( i = 1 ; i <= Mobile_Structure.length ; i ++ ) {
    free( Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Mobile_Structure.Residue ) ;

/************/

  printf( "Finished\n\n" ) ;

  return( 0 ) ;

} /* end main */
