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

  int	i , j , k ;

  /* Command line options */

  char		*input_file_name ;
  char		*output_file_name ;
  float		distance ;
  char		*constraints ;
  int		number_of_constraints ;
  char		*constraint[NUMBER_OF_CONSTRAINTS] ;

  char		*default_distance ;
  char		*default_constraints ;

  char		*gf_type ;
  char		gf_one_letter_code ;

  /* File stuff */

  char		*static_file_name ;
  char		*mobile_file_name ;

  FILE		*ftdock_in_file ;
  FILE		*ftdock_out_file ;
  char		line_buffer[100] ;
  char		saved_header_line[SAVED_HEADER_LINES][100] ;
  int		saved_header_lines ;

  /* ftdock.dat file values */

  int		complex_id , prev_complex_id , new_complex_id , SCscore ;
  float		RPscore ;

  /* Angles stuff */

  int		z_twist , theta , phi ;

  /* Structures */

  struct Structure	Static_Structure , Mobile_Structure ;
  struct Structure	Origin_Static_Structure , Origin_Mobile_Structure ;
  struct Structure	Sparse_Origin_Static_Structure , Sparse_Origin_Mobile_Structure ;
  struct Structure	Rotated_Sparse_Origin_Mobile_Structure , Translated_and_Rotated_Sparse_Origin_Mobile_Structure ;

  int		useful ;
  char		*half_constraint[NUMBER_OF_CONSTRAINTS * 2] ;
  char		*residue_number_plus_iCode ;
  int		static_residue , static_atom ;
  int		mobile_residue , mobile_atom ;
  int		prox ;

  /* Co-ordinates */

  int		x , y , z ;

  /* Grid stuff */

  int		grid_size ;
  float		grid_span ;
  float		one_span ;

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

  printf( "Starting Filter program\n" ) ;

/************/

  /* Memory allocation */

  if( ( ( input_file_name   = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( output_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( constraints       = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( gf_type           = ( char * ) malloc (  10 * sizeof( char ) ) ) == NULL ) ||
      ( ( static_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( mobile_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ) {
    GENERAL_MEMORY_PROBLEM 
  }

  for( i = 0 ; i < NUMBER_OF_CONSTRAINTS ; i ++ ) {
    if( ( constraint[i]  = ( char * ) malloc ( 50 * sizeof( char ) ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }
  }
  for( i = 0 ; i < NUMBER_OF_CONSTRAINTS * 2 ; i ++ ) {
    if( ( half_constraint[i]  = ( char * ) malloc ( 50 * sizeof( char ) ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }
  }

  if( ( residue_number_plus_iCode = ( char * ) malloc ( 50 * sizeof( char ) ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

/************/

  /* Command Line defaults */

  strcpy( input_file_name , "ftdock_rpscored.dat" ) ;
  strcpy( output_file_name , "ftdock_filtered.dat" ) ;
  distance = 4.5 ;
  strcpy( constraints , "" ) ;

  default_distance = "(default)" ;
  default_constraints = "error" ;
  number_of_constraints = 0 ;

  strcpy( gf_type , "G_DATA" ) ;
  gf_one_letter_code = 'g' ;

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
        if( strcmp( argv[i] , "-distance" ) == 0 ) {
          i ++ ;
          if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
          }
          sscanf( argv[i] , "%f" , &distance ) ;
          default_distance = "(user defined)" ;
        } else {
          if( strcmp( argv[i] , "-constraints" ) == 0 ) {
            i ++ ;
            number_of_constraints = -1 ;
            if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
              printf( "Bad command line\n" ) ;
              exit( EXIT_FAILURE ) ;
            }
            strcat( constraints , argv[i] ) ;
            strcpy( constraint[++number_of_constraints] , argv[i] ) ;
            while( ( (i+1) < argc ) && ( strncmp( argv[i+1] , "-" , 1 ) != 0 ) ) {
              if( number_of_constraints == ( NUMBER_OF_CONSTRAINTS - 1 ) ) {
                printf( "Too many constraints\nDying\n" ) ;
                exit( EXIT_FAILURE ) ;
              }
              i ++ ;
              strcat( constraints , " " ) ;
              strcat( constraints , argv[i] ) ;
              strcat( constraint[++number_of_constraints] , argv[i] ) ;
            }
            default_constraints = "defined" ;
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

  }

  if( strcmp( default_constraints , "error" ) == 0 ) {
    printf( "No constraints provided - pointless\nDying\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

  if( strcmp( input_file_name , output_file_name ) == 0 ) {
    printf( "Am not happy with input and output files being the same name\nDying\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

/************/

  /* Read data file - headers */

  if( ( ftdock_in_file = fopen( input_file_name , "r" ) ) == NULL ) {
    printf( "Could not open %s for reading.\nDying\n\n" , input_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }

  saved_header_lines = 0 ;

  while( fgets( line_buffer , 99 , ftdock_in_file ) ) {

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

  /* Make the structures sparse - no use wasting time with junk */

  /* Parse the constraints */

  for( i = 0 ; i <= number_of_constraints ; i ++ ) {
    strcpy( half_constraint[i+number_of_constraints+1] , strchr( constraint[i] , ':' ) ) ;
    half_constraint[i+number_of_constraints+1] ++ ;
    strcpy( strchr( constraint[i] , ':' ) , "\0" ) ;
    strcpy( half_constraint[i] , constraint[i] ) ;
  }

/************/

  /* Static Structure */

  if( ( Sparse_Origin_Static_Structure.Residue = ( struct Amino_Acid * ) malloc ( ( Origin_Static_Structure.length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  strcpy( Sparse_Origin_Static_Structure.ident , "Complex" ) ;
  j = 0 ;
  for( i = 1 ; i <= Origin_Static_Structure.length ; i ++ ) {
    useful = 0 ;
    for( k = 0 ; k <= 1 + number_of_constraints * 2 ; k ++ ) {
      if( strpbrk( half_constraint[k] , Origin_Static_Structure.Residue[i].chainID ) != NULL ) {
        if( strpbrk( half_constraint[k] , "0123456789" ) != NULL ) {
          strcpy( residue_number_plus_iCode , strpbrk( half_constraint[k] , "0123456789" ) ) ;
          if( strncmp( strpbrk( Origin_Static_Structure.Residue[i].res_seq_plus_iCode , "0123456789" ) ,  strpbrk( residue_number_plus_iCode , "0123456789" ) , (int)strlen( residue_number_plus_iCode ) ) == 0 ) { useful = 1 ; }
        } else {
          useful = 1 ;
        }
      }
    }
    if( useful != 0 ) {
      j ++ ;
      if( ( Sparse_Origin_Static_Structure.Residue[j].Atom = ( struct Atom * ) malloc ( ( Origin_Static_Structure.Residue[i].size +1 ) * sizeof_Atom ) ) == NULL ) {
        GENERAL_MEMORY_PROBLEM
      }
      strcpy( Sparse_Origin_Static_Structure.Residue[j].res_name           , Origin_Static_Structure.Residue[i].res_name ) ;
      strcpy( Sparse_Origin_Static_Structure.Residue[j].chainID            , Origin_Static_Structure.Residue[i].chainID ) ;
      strcpy( Sparse_Origin_Static_Structure.Residue[j].res_seq_plus_iCode , Origin_Static_Structure.Residue[i].res_seq_plus_iCode ) ;
      strcpy( Sparse_Origin_Static_Structure.Residue[j].olc                , Origin_Static_Structure.Residue[i].olc ) ;
      Sparse_Origin_Static_Structure.Residue[j].nc                 = Origin_Static_Structure.Residue[i].nc   ;
      Sparse_Origin_Static_Structure.Residue[j].size               = Origin_Static_Structure.Residue[i].size ;
      for( k = 1 ; k <= Origin_Static_Structure.Residue[i].size ; k ++ ) {
        Sparse_Origin_Static_Structure.Residue[j].Atom[k] = Origin_Static_Structure.Residue[i].Atom[k] ;
      }
    }
  }
  Sparse_Origin_Static_Structure.length = j ;

  if( ( Sparse_Origin_Static_Structure.Residue = ( struct Amino_Acid * ) realloc ( Sparse_Origin_Static_Structure.Residue , ( Sparse_Origin_Static_Structure.length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  
/************/

  /* Mobile Structure */

  if( ( Sparse_Origin_Mobile_Structure.Residue = ( struct Amino_Acid * ) malloc ( ( Origin_Mobile_Structure.length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  strcpy( Sparse_Origin_Mobile_Structure.ident , "Complex" ) ;
  j = 0 ;
  for( i = 1 ; i <= Origin_Mobile_Structure.length ; i ++ ) {
    useful = 0 ;
    for( k = 0 ; k <= 1 + number_of_constraints * 2 ; k ++ ) {
      if( strpbrk( half_constraint[k] , Origin_Mobile_Structure.Residue[i].chainID ) != NULL ) {
        if( strpbrk( half_constraint[k] , "0123456789" ) != NULL ) {
          strcpy( residue_number_plus_iCode , strpbrk( half_constraint[k] , "0123456789" ) ) ;
          if( strncmp( strpbrk( Origin_Mobile_Structure.Residue[i].res_seq_plus_iCode , "0123456789" ) ,  strpbrk( residue_number_plus_iCode , "0123456789" ) , (int)strlen( residue_number_plus_iCode ) ) == 0 ) { useful = 1 ; }
        } else {
          useful = 1 ;
        }
      }
    }
    if( useful != 0 ) {
      j ++ ;
      if( ( Sparse_Origin_Mobile_Structure.Residue[j].Atom = ( struct Atom * ) malloc ( ( Origin_Mobile_Structure.Residue[i].size +1 ) * sizeof_Atom ) ) == NULL ) {
        GENERAL_MEMORY_PROBLEM
      }
      strcpy( Sparse_Origin_Mobile_Structure.Residue[j].res_name           , Origin_Mobile_Structure.Residue[i].res_name ) ;
      strcpy( Sparse_Origin_Mobile_Structure.Residue[j].chainID            , Origin_Mobile_Structure.Residue[i].chainID ) ;
      strcpy( Sparse_Origin_Mobile_Structure.Residue[j].res_seq_plus_iCode , Origin_Mobile_Structure.Residue[i].res_seq_plus_iCode ) ;
      strcpy( Sparse_Origin_Mobile_Structure.Residue[j].olc                , Origin_Mobile_Structure.Residue[i].olc ) ;
      Sparse_Origin_Mobile_Structure.Residue[j].nc                 = Origin_Mobile_Structure.Residue[i].nc   ;
      Sparse_Origin_Mobile_Structure.Residue[j].size               = Origin_Mobile_Structure.Residue[i].size ;
      for( k = 1 ; k <= Origin_Mobile_Structure.Residue[i].size ; k ++ ) {
        Sparse_Origin_Mobile_Structure.Residue[j].Atom[k] = Origin_Mobile_Structure.Residue[i].Atom[k] ;
      }
    }
  }
  Sparse_Origin_Mobile_Structure.length = j ;

  if( ( Sparse_Origin_Mobile_Structure.Residue = ( struct Amino_Acid * ) realloc ( Sparse_Origin_Mobile_Structure.Residue , ( Sparse_Origin_Mobile_Structure.length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

/************/

  /* Free some memory */

  for( i = 1 ; i <= Origin_Static_Structure.length ; i ++ ) {
    free( Origin_Static_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Static_Structure.Residue ) ;
  for( i = 1 ; i <= Origin_Mobile_Structure.length ; i ++ ) {
    free( Origin_Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Mobile_Structure.Residue ) ;

/************/

  /* Main program loop - some initialisations */

  printf( "Starting main loop\n" ) ;

  /* Open output ftdock.dat */

  if( ( ftdock_out_file = fopen( output_file_name , "w" ) ) == NULL ) {
    printf( "Could not open %s for writing.\nDying\n\n" , output_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }

  setvbuf( ftdock_out_file , (char *)NULL , _IONBF , 0 ) ;

  for( i = 1 ; i <= saved_header_lines ; i ++ ) {
    fprintf( ftdock_out_file , saved_header_line[i] ) ;
  }

  fprintf( ftdock_out_file, "******\n\n" ) ;

  fprintf( ftdock_out_file, "Filter\n" ) ;

  /* Write out command line options */

  fprintf( ftdock_out_file, "\nCommand line controllable values\n" ) ;
  fprintf( ftdock_out_file, "Constraints                        :: %s ( %d )\n" , constraints , number_of_constraints + 1 ) ;
  fprintf( ftdock_out_file, "Distance                           :: %9.2f   %s\n" , distance , default_distance ) ;

  fprintf( ftdock_out_file, "\nData\n" ) ;
  fprintf( ftdock_out_file , "Type       ID    prvID    SCscore        RPscore         Coordinates            Angles\n\n" ) ;

  new_complex_id = 0 ;

/************/

  /* Main program loop */

  /* No buffering at all */
  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "Filtering ( one dot / complex ) " ) ;

  while( fgets( line_buffer, 99, ftdock_in_file ) ) {

    if( strncmp( line_buffer , gf_type , 6 ) == 0 ) {

      sscanf( line_buffer + 1 , "_DATA %d %d %d %f  %d %d %d  %d %d %d" , &complex_id , &prev_complex_id , &SCscore , &RPscore , &x , &y , &z , &z_twist , &theta , &phi ) ;

      printf( "." ) ;

      /* Build sparse molecule */

      /* Rotate Mobile Structure */
      Rotated_Sparse_Origin_Mobile_Structure = rotate_structure( Sparse_Origin_Mobile_Structure , z_twist , theta , phi ) ;

      /* Translate Mobile Structure */
      Translated_and_Rotated_Sparse_Origin_Mobile_Structure = translate_structure( Rotated_Sparse_Origin_Mobile_Structure , (float)x * one_span , (float)y * one_span , (float)z * one_span ) ;

      /* Calculate all atom--atom distances */

      prox = 0 ;

      for( static_residue = 1 ; static_residue <= Sparse_Origin_Static_Structure.length ; static_residue ++ ) {
        for( mobile_residue = 1 ; mobile_residue <= Sparse_Origin_Mobile_Structure.length ; mobile_residue ++ ) {

          for( static_atom = 1 ; static_atom <= Sparse_Origin_Static_Structure.Residue[static_residue].size ; static_atom ++ ) {
            for( mobile_atom = 1 ; mobile_atom <= Sparse_Origin_Mobile_Structure.Residue[mobile_residue].size ; mobile_atom ++ ) {

              if( pythagoras( Sparse_Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[1] , Sparse_Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[2] , Sparse_Origin_Static_Structure.Residue[static_residue].Atom[static_atom].coord[3] , Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[1] , Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[2] , Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[mobile_residue].Atom[mobile_atom].coord[3] ) < distance ) {

                for( k = 0 ; k <= number_of_constraints ; k ++ ) {
                  useful = 0 ;
                  if( strpbrk( half_constraint[k] , Sparse_Origin_Static_Structure.Residue[static_residue].chainID ) != NULL ) {
                    if( strpbrk( half_constraint[k] , "0123456789" ) != NULL ) {
                      strcpy( residue_number_plus_iCode , strpbrk( half_constraint[k] , "0123456789" ) ) ;
                      if( strncmp( strpbrk( Sparse_Origin_Static_Structure.Residue[static_residue].res_seq_plus_iCode , "0123456789" ) ,  strpbrk( residue_number_plus_iCode , "0123456789" ) , (int)strlen( residue_number_plus_iCode ) ) == 0 ) { useful = 1 ; }
                    } else {
                      useful = 1 ;
                    }
                  }
                  if( useful == 1 ) {
                    if( strpbrk( half_constraint[k+number_of_constraints+1] , Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[mobile_residue].chainID ) != NULL ) {
                      if( strpbrk( half_constraint[k+number_of_constraints+1] , "0123456789" ) != NULL ) {
                        strcpy( residue_number_plus_iCode , strpbrk( half_constraint[k+number_of_constraints+1] , "0123456789" ) ) ;
                        if( strncmp( strpbrk( Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[mobile_residue].res_seq_plus_iCode , "0123456789" ) ,  strpbrk( residue_number_plus_iCode , "0123456789" ) , (int)strlen( residue_number_plus_iCode ) ) == 0 ) { prox = 1 ; }
                      } else {
                        prox = 1 ;
                      }
                    }
                    
                  }
                }
                for( k = 0 ; k <= number_of_constraints ; k ++ ) {
                  useful = 0 ;
                  if( strpbrk( half_constraint[k] , Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[mobile_residue].chainID ) != NULL ) {
                    if( strpbrk( half_constraint[k] , "0123456789" ) != NULL ) {
                      strcpy( residue_number_plus_iCode , strpbrk( half_constraint[k] , "0123456789" ) ) ;
                      if( strncmp( strpbrk( Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[mobile_residue].res_seq_plus_iCode , "0123456789" ) ,  strpbrk( residue_number_plus_iCode , "0123456789" ) , (int)strlen( residue_number_plus_iCode ) ) == 0 ) { useful = 1 ; }
                    } else {
                      useful = 1 ;
                    }
                  }
                  if( useful == 1 ) {
                    if( strpbrk( half_constraint[k+number_of_constraints+1] , Sparse_Origin_Static_Structure.Residue[static_residue].chainID ) != NULL ) {
                      if( strpbrk( half_constraint[k+number_of_constraints+1] , "0123456789" ) != NULL ) {
                        strcpy( residue_number_plus_iCode , strpbrk( half_constraint[k+number_of_constraints+1] , "0123456789" ) ) ;
                        if( strncmp( strpbrk( Sparse_Origin_Static_Structure.Residue[static_residue].res_seq_plus_iCode , "0123456789" ) ,  strpbrk( residue_number_plus_iCode , "0123456789" ) , (int)strlen( residue_number_plus_iCode ) ) == 0 ) { prox = 1 ; }
                      } else {
                        prox = 1 ;
                      }
                    }
                  
                  }
                }

              }
            }

          }

        }

      }

      if( prox == 1 ) {

        new_complex_id ++ ;

        fprintf( ftdock_out_file, "%6s %6d   %6d    %7d       %8.3f      %4d %4d %4d      %4d%4d%4d\n" , gf_type , new_complex_id , complex_id , SCscore , RPscore , x , y , z , z_twist , theta , phi ) ;

      }

      /* Free some memory */
      for( j = 1 ; j <= Rotated_Sparse_Origin_Mobile_Structure.length ; j ++ ) {
        free( Rotated_Sparse_Origin_Mobile_Structure.Residue[j].Atom ) ;
        free( Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue[j].Atom ) ;
      }
      free( Rotated_Sparse_Origin_Mobile_Structure.Residue ) ;
      free( Translated_and_Rotated_Sparse_Origin_Mobile_Structure.Residue ) ;

    }
    
  }

  /* Finished main loop */

  fclose( ftdock_in_file ) ;
  fclose( ftdock_out_file ) ;

  printf( "\n\n" ) ;

/************/

  /* Free the memory */

  for( i = 1 ; i <= Sparse_Origin_Static_Structure.length ; i ++ ) {
    free( Sparse_Origin_Static_Structure.Residue[i].Atom ) ;
  }
  free( Sparse_Origin_Static_Structure.Residue ) ;

  for( i = 1 ; i <= Sparse_Origin_Mobile_Structure.length ; i ++ ) {
    free( Sparse_Origin_Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Sparse_Origin_Mobile_Structure.Residue ) ;

/************/

  printf( "Finished\n\n" ) ;

  return( 0 ) ;

} /* end main */
