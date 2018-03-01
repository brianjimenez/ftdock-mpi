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

  int		i , j ;

  /* Command line options */

  char		*input_file_name ;
  char		*output_file_name ;
  char		*matrix_file_name ;
  char		*gf_type ;
  char		gf_one_letter_code ;

  char		*default_matrix_file_name ;

  /* File stuff */

  char		*static_file_name ;
  char		*mobile_file_name ;

  FILE		*ftdock_in_file ;
  FILE		*ftdock_out_file ;
  FILE		*matrix_file ;
  char		line_buffer[200] ;
  char		saved_header_line[SAVED_HEADER_LINES][100] ;
  int		saved_header_lines ;

  /* ftdock_*.dat file values */

  int		complex_id , prev_complex_id , SCscore ;
  float		RPscore ;

  /* Angles stuff */

  int		z_twist , theta , phi ;

  /* Structures */

  struct Structure	Static_Structure , Mobile_Structure ;
  struct Structure	Origin_Static_Structure , Origin_Mobile_Structure ;
  struct Structure	Rotated_Origin_Mobile_Structure , Translated_and_Rotated_Origin_Mobile_Structure ;

  /* Co-ordinates */

  int		x , y , z ;

  /* Grid stuff */

  int		grid_size ;
  float		grid_span ;
  float		one_span ;

  int		static_residue , static_atom ;
  int		mobile_residue , mobile_atom ;
  int		prox ;

  struct Matrix	Matrix ;

  /* Scores */

  int		interesting_scores ;
  struct Score	*Scores ;

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

  printf( "Starting RPScore residue level pair potential score program\n" ) ;

/************/

  /* Memory allocation */

  if( ( ( input_file_name   = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( output_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( matrix_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( gf_type           = ( char * ) malloc (  10 * sizeof( char ) ) ) == NULL ) ||
      ( ( static_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( mobile_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ) {
    GENERAL_MEMORY_PROBLEM 
  }

/************/

  /* Command Line defaults */

  strcpy( input_file_name , "ftdock_global.dat" ) ;
  strcpy( output_file_name , "ftdock_rpscored.dat" ) ;
  strcpy( matrix_file_name , "best.matrix" ) ;

  strcpy( gf_type , "G_DATA" ) ;
  gf_one_letter_code = 'g' ;

  default_matrix_file_name = "(default)" ;

  /* Command Line parse */

  for( i = 1 ; i < argc ; i ++ ) {


    if( strcmp( argv[i] , "-in" ) == 0 ) {
      i ++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
        printf( "Bad command line\n" ) ;
        exit( EXIT_FAILURE ) ;
      }
      strcpy( input_file_name , argv[i] ) ;
    } else {
      if( strcmp( argv[i] , "-out" ) == 0 ) {
        i ++ ;
        if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
          printf( "Bad command line\n" ) ;
          exit( EXIT_FAILURE ) ;
        }
        strcpy( output_file_name , argv[i] ) ;
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
          if( strcmp( argv[i] , "-fine" ) == 0 ) {
            strcpy( gf_type , "F_DATA" ) ;
            gf_one_letter_code = 'f' ;
          } else {
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
          }
        }
      }
    }

  }

  if( strcmp( input_file_name , output_file_name ) == 0 ) {
    printf( "Am not happy with input and output files being the same name\nDying\n" ) ;
    exit( EXIT_FAILURE ) ;
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

  /* Read ftdock file - headers */

  if( ( ftdock_in_file = fopen( input_file_name , "r" ) ) == NULL ) {
    printf( "Could not open %s for reading.\nDying\n\n" , input_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }

  saved_header_lines = 0 ;

  while( fgets( line_buffer , 199 , ftdock_in_file ) ) {

    if( strncmp( line_buffer , "Static molecule" , 15 ) == 0 ) sscanf( line_buffer , "Static molecule                :: %s" , static_file_name ) ;
    if( strncmp( line_buffer , "Mobile molecule" , 15 ) == 0 ) sscanf( line_buffer , "Mobile molecule                :: %s" , mobile_file_name ) ;
    if( gf_one_letter_code == 'g' ) {
      if( strncmp( line_buffer , "Global grid size" , 16 ) == 0 ) sscanf( line_buffer , "Global grid size               :: %d" , &grid_size ) ;
    } else {
      if( strncmp( line_buffer , "Refinement grid size" , 20 ) == 0 ) sscanf( line_buffer , "Refinement grid size               :: %d" , &grid_size ) ;
    }

    if( strncmp( line_buffer , "Data" , 4 ) == 0 ) break ;

    strncpy( saved_header_line[++saved_header_lines] , line_buffer , 99 ) ;

  }

/************/

  /* Read in Structures from pdb files */
  Static_Structure = read_pdb_to_structure( static_file_name ) ;
  Mobile_Structure = read_pdb_to_structure( mobile_file_name ) ;

/************/

  /* Store new structures centered on Origin */

  Origin_Static_Structure = translate_structure_onto_origin( Static_Structure ) ;
  Origin_Mobile_Structure = translate_structure_onto_origin( Mobile_Structure ) ;

  /* Free some memory */

  for( i = 1 ; i <= Static_Structure.length ; i ++ ) {
    free( Static_Structure.Residue[i].Atom ) ;
  }
  free( Static_Structure.Residue ) ;
  for( i = 1 ; i <= Mobile_Structure.length ; i ++ ) {
    free( Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Mobile_Structure.Residue ) ;

/************/

  /* Calculate Grid stuff */

  grid_span = total_span_of_structures( Origin_Static_Structure , Origin_Mobile_Structure ) ;

  one_span = grid_span / (float)grid_size ;

/************/

  /* Main program loop - some initialisations */

  /* Open output data file */

  if( ( ftdock_out_file = fopen( output_file_name , "w" ) ) == NULL ) {
    printf( "Could not open %s for writing.\nDying\n\n" , output_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }
  /* Line buffering only - so user can see whats going on a bit better */
  setvbuf( ftdock_out_file , (char *)NULL , _IOLBF , 0 ) ;

  for( i = 1 ; i <= saved_header_lines ; i ++ ) {
    fprintf( ftdock_out_file , saved_header_line[i] ) ;
  }

  fprintf( ftdock_out_file, "******\n\n" ) ;

  fprintf( ftdock_out_file, "Residue level Pair Potential Scoring\n" ) ;

  /* Write out command line options */

  fprintf( ftdock_out_file, "\nCommand line controllable values\n" ) ;
  fprintf( ftdock_out_file, "Matrix                             :: %s    %s\n" , matrix_file_name , default_matrix_file_name ) ;

  fprintf( ftdock_out_file, "\nData\n" ) ;
  fprintf( ftdock_out_file , "Type       ID    prvID    SCscore        RPscore         Coordinates            Angles\n\n" ) ;

  interesting_scores = 0 ;

  fclose( ftdock_out_file ) ;

  if( ( ftdock_out_file = fopen( "rpscore_scratch_scores" , "w" ) ) == NULL ) {
    printf( "Could not open scratch file for writing.\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }
  /* Line buffering only - so user can see whats going on a bit better */
  setvbuf( ftdock_out_file , (char *)NULL , _IOLBF , 0 ) ;

/************/

  /* Main program loop */

  /* No buffering at all */
  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "Scoring using Residue level Pair Potentials ( one dot / complex ) " ) ;

  while( fgets( line_buffer, 99, ftdock_in_file ) ) {

    if( strncmp( line_buffer , gf_type , 6 ) == 0 ) {

      sscanf( line_buffer + 1 , "_DATA %d %d %d %f  %d %d %d  %d %d %d" , &complex_id , &prev_complex_id , &SCscore , &RPscore , &x , &y , &z , &z_twist , &theta , &phi ) ;

      RPscore = 0 ;

      printf( "." ) ;

      /* Rotate Mobile Structure */
      Rotated_Origin_Mobile_Structure = rotate_structure( Origin_Mobile_Structure , z_twist , theta , phi ) ;

      /* Translate Mobile Structure */
      Translated_and_Rotated_Origin_Mobile_Structure = translate_structure( Rotated_Origin_Mobile_Structure , (float)x * one_span , (float)y * one_span , (float)z * one_span ) ;

      for( static_residue = 1 ; static_residue <= Origin_Static_Structure.length ; static_residue ++ ) {

        if( Origin_Static_Structure.Residue[static_residue].nc > 20 ) break ;

        for( mobile_residue = 1 ; mobile_residue <= Origin_Mobile_Structure.length ; mobile_residue ++ ) {

          if( Origin_Mobile_Structure.Residue[mobile_residue].nc > 20 ) break ;
          prox = 0 ;

          for( static_atom = 1 ; static_atom <= Origin_Static_Structure.Residue[static_residue].size ; static_atom ++ ) {
            for( mobile_atom = 1 ; mobile_atom <= Origin_Mobile_Structure.Residue[mobile_residue].size ; mobile_atom ++ ) {

              if( ( ( Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[1] - Translated_and_Rotated_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[1] ) * ( Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[1] - Translated_and_Rotated_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[1] ) ) + ( ( Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[2] - Translated_and_Rotated_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[2] ) * ( Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[2] - Translated_and_Rotated_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[2] ) ) + ( ( Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[3] - Translated_and_Rotated_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[3] ) * ( Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[3] - Translated_and_Rotated_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[3] ) ) < ( Matrix.distance * Matrix.distance ) ) prox = 1 ;

            }
          }
          if( prox == 1 ) RPscore += Matrix.score[Origin_Static_Structure.Residue[static_residue].nc][Origin_Mobile_Structure.Residue[mobile_residue].nc] ;

        }
      }

      interesting_scores ++ ;

      fprintf( ftdock_out_file, "%6s %5d   %5d    %7d       %8.3f      %4d %4d %4d      %4d%4d%4d\n" , gf_type , complex_id , complex_id , SCscore , RPscore , x , y , z , z_twist , theta , phi ) ;

      /* Free some memory */
      for( j = 1 ; j <= Origin_Mobile_Structure.length ; j ++ ) {
        free( Rotated_Origin_Mobile_Structure.Residue[j].Atom ) ;
        free( Translated_and_Rotated_Origin_Mobile_Structure.Residue[j].Atom ) ;
      }
      free( Rotated_Origin_Mobile_Structure.Residue ) ;
      free( Translated_and_Rotated_Origin_Mobile_Structure.Residue ) ;

    }

  }

  /* Finished main loop */

  fclose( ftdock_in_file ) ;
  fclose( ftdock_out_file ) ;

  printf( "\n" ) ;

/************/

  /* Free the memory */

  for( i = 1 ; i <= Origin_Static_Structure.length ; i ++ ) {
    free( Origin_Static_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Static_Structure.Residue ) ;

  for( i = 1 ; i <= Origin_Mobile_Structure.length ; i ++ ) {
    free( Origin_Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Mobile_Structure.Residue ) ;

/************/

  /* Read in all the scores */

  printf( "\nReloading all the scores\n" ) ;

  if( ( ftdock_in_file = fopen( "rpscore_scratch_scores" , "r" ) ) == NULL ) {
    printf( "Could not open scratch file for readng.\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

  if( ( Scores = ( struct Score * ) malloc ( interesting_scores * sizeof( struct Score ) ) ) == NULL ) {
    printf( "Not enough memory left for storing scores\nProbably keeping too many per rotation\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

  interesting_scores = 0 ;

  while( fgets( line_buffer , 99 , ftdock_in_file ) ) {

    sscanf( line_buffer + 1 , "_DATA %d %d %d %f  %d %d %d  %d %d %d" , &complex_id , &prev_complex_id , &SCscore , &RPscore ,
                                                                     &x , &y , &z , &z_twist , &theta , &phi ) ;

    Scores[interesting_scores].score    = SCscore ;
    Scores[interesting_scores].coord[1] = x ;
    Scores[interesting_scores].coord[2] = y ;
    Scores[interesting_scores].coord[3] = z ;
    Scores[interesting_scores].angle[1] = z_twist ;
    Scores[interesting_scores].angle[2] = theta ;
    Scores[interesting_scores].angle[3] = phi ;
    Scores[interesting_scores].rpscore  = RPscore ;
    Scores[interesting_scores].extra    = complex_id ;

    interesting_scores ++ ;

  }

  fclose( ftdock_in_file ) ;

  interesting_scores -- ;

  qsort_rpscores( Scores , 0 , interesting_scores ) ;

/************/

  if( ( ftdock_out_file = fopen( output_file_name , "a" ) ) == NULL ) {
    printf( "Could not open %s for writing.\nDying\n\n" , output_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }
  /* Line buffering only - so user can see whats going on a bit better */
  setvbuf( ftdock_out_file , (char *)NULL , _IOLBF , 0 ) ;

  for( i = 0 ; i <= interesting_scores ; i ++ ) {

    fprintf( ftdock_out_file, "%6s %6d   %6d    %7d       %8.3f      %4d %4d %4d      %4d%4d%4d\n" , gf_type ,
               i + 1 , Scores[i].extra , Scores[i].score , Scores[i].rpscore ,
               Scores[i].coord[1] , Scores[i].coord[2] , Scores[i].coord[3] ,
               Scores[i].angle[1] , Scores[i].angle[2] , Scores[i].angle[3] ) ;

  }

  fclose( ftdock_out_file ) ;
    
/************/

  printf( "Finished\n\n" ) ;

  return( 0 ) ;

} /* end main */
