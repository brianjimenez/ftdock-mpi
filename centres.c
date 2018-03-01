/*
Copyright (C) 1997-2000 Imperial Cancer Research Technology
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

  int	i , j ;

  int	residue , atom , total_atoms ;
  int	length , atom_serial_number , residue_number ;

  /* Command line options */

  char		*input_file_name ;
  char		*gf_type ;
  char		gf_one_letter_code ;

  /* File stuff */

  FILE		*ftdock_file ;
  char		line_buffer[100] ;

  char		*static_file_name ;
  char		*mobile_file_name ;

  /* Structures */

  struct Structure	Static_Structure , Mobile_Structure ;
  struct Structure	Centres_Structure , Origin_Mobile_Structure ;
  struct Structure	Rotated_at_Origin_Mobile_Structure , Translated_and_Rotated_at_Origin_Mobile_Structure ;

  /* Co-ordinates and Angles */

  int	x , y , z , z_twist , theta, phi ;

  float		average_x , average_y , average_z ;

  /* Grid stuff */

  int	grid_size ;
  float		grid_span , one_span ;

  /* ftdock.dat file values */

  int	complex_id , prev_complex_id , SCscore ;
  float		RPscore ;

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

  printf( "Starting Centres central positions of mobile molecules program\n" ) ;

/************/

  /* Memory allocation */

  if( ( ( input_file_name   = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( gf_type           = ( char * ) malloc (  10 * sizeof( char ) ) ) == NULL ) ||
      ( ( static_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( mobile_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ) {
    GENERAL_MEMORY_PROBLEM
  }

/************/

  /* Command Line defaults */

  strcpy( input_file_name , "ftdock_global.dat" ) ;
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
      if( strcmp( argv[i] , "-fine" ) == 0 ) {
        strcpy( gf_type , "F_DATA" ) ;
        gf_one_letter_code = 'f' ;
      } else {
        printf( "Bad command line\n" ) ;
        exit( EXIT_FAILURE ) ;
      }
    }

  }

/************/

  /* Read neat data file - headers */

  if( ( ftdock_file = fopen( input_file_name , "r" ) ) == NULL ) {
    printf( "Could not open %s for reading.\nDying\n\n" , input_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }

  while( fgets( line_buffer, 99, ftdock_file ) ) {

    if( strncmp( line_buffer , "Static molecule" , 15 ) == 0 ) sscanf( line_buffer , "Static molecule :: %s" , static_file_name ) ;

    if( strncmp( line_buffer , "Mobile molecule" , 15 ) == 0 ) sscanf( line_buffer , "Mobile molecule :: %s" , mobile_file_name ) ;

    if( gf_one_letter_code == 'g' ) {

      if( strncmp( line_buffer , "Global grid size" , 16 ) == 0 ) sscanf( line_buffer , "Global grid size :: %d" , &grid_size ) ;

    } else {

      if( strncmp( line_buffer , "Refinement grid size" , 20 ) == 0 ) sscanf( line_buffer , "Refinement grid size :: %d" , &grid_size ) ;

    }

    if( strncmp( line_buffer , gf_type , 6 ) == 0 ) break ;

  }

  fclose( ftdock_file ) ;

/************/

  /* Read in Structures from pdb files */
  Static_Structure = read_pdb_to_structure( static_file_name ) ;
  Mobile_Structure = read_pdb_to_structure( mobile_file_name ) ;

/************/

  /* Store new structures centered on Origin */

  Centres_Structure = translate_structure_onto_origin( Static_Structure ) ;
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

  grid_span = total_span_of_structures( Centres_Structure , Origin_Mobile_Structure ) ;

  one_span = grid_span / (float)grid_size ;

/************/

  /* Reopen ftdock.dat */

  if( ( ftdock_file = fopen( input_file_name , "r" ) ) == NULL ) {
    printf( "Could not open %s for reading.\nDying\n\n" , input_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }

/************/

  /* Loop */

  length = Centres_Structure.length ;
  atom_serial_number = Centres_Structure.Residue[length].Atom[1].serial + 100 ;
  residue_number = 0 ;

  while( fgets( line_buffer, 99, ftdock_file ) ) {

    if( strncmp( line_buffer , gf_type , 6 ) == 0 ) {

      sscanf( line_buffer + 1 , "_DATA %d %d %d %f  %d %d %d  %d %d %d" , &complex_id , &prev_complex_id , &SCscore , &RPscore , &x , &y , &z , &z_twist , &theta , &phi ) ;

      /* Rotate Mobile Structure */
      Rotated_at_Origin_Mobile_Structure = rotate_structure( Origin_Mobile_Structure , z_twist , theta , phi ) ;
         
      /* Translate Mobile Structure */
      Translated_and_Rotated_at_Origin_Mobile_Structure = translate_structure( Rotated_at_Origin_Mobile_Structure , (float)x * one_span , (float)y * one_span , (float)z * one_span ) ;
         
      /* Find current centre of mobile structure */

      total_atoms = 0 ;

      average_x = 0 ;
      average_y = 0 ;
      average_z = 0 ;

      for( residue = 1 ; residue <= Translated_and_Rotated_at_Origin_Mobile_Structure.length ; residue ++ ) {

        for( atom = 1 ; atom <= Translated_and_Rotated_at_Origin_Mobile_Structure.Residue[residue].size ; atom ++ ) {

          total_atoms ++ ;

          average_x += Translated_and_Rotated_at_Origin_Mobile_Structure.Residue[residue].Atom[atom].coord[1] ;
          average_y += Translated_and_Rotated_at_Origin_Mobile_Structure.Residue[residue].Atom[atom].coord[2] ;
          average_z += Translated_and_Rotated_at_Origin_Mobile_Structure.Residue[residue].Atom[atom].coord[3] ;

        }

      }

      average_x /= (float)total_atoms ;
      average_y /= (float)total_atoms ;
      average_z /= (float)total_atoms ;

      /* Add centre */

      length ++ ;
      atom_serial_number ++ ;
      residue_number ++ ;

      if( ( Centres_Structure.Residue = (struct Amino_Acid * ) realloc ( Centres_Structure.Residue, ( length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
        GENERAL_MEMORY_PROBLEM
      }

      strcpy( Centres_Structure.Residue[length].chainID            , "X" ) ;
      sprintf( Centres_Structure.Residue[length].res_seq_plus_iCode , "%4d " , residue_number ) ;
      strcpy( Centres_Structure.Residue[length].olc                , "x" ) ;
      Centres_Structure.Residue[length].nc                         = 0 ;
      Centres_Structure.Residue[length].size                       = 1 ;
      strcpy( Centres_Structure.Residue[length].res_name           , "HOH" ) ;

      if( ( Centres_Structure.Residue[length].Atom = ( struct Atom * ) malloc ( 2 * sizeof_Atom ) ) == NULL ) {
        GENERAL_MEMORY_PROBLEM
      }
      Centres_Structure.Residue[length].Atom[1].serial = atom_serial_number ;
      strcpy( Centres_Structure.Residue[length].Atom[1].atom_name, "  O " ) ;
      Centres_Structure.Residue[length].Atom[1].coord[1] = average_x ;
      Centres_Structure.Residue[length].Atom[1].coord[2] = average_y ;
      Centres_Structure.Residue[length].Atom[1].coord[3] = average_z ;
      Centres_Structure.Residue[length].Atom[1].occupancy = 0.0 ;
      Centres_Structure.Residue[length].Atom[1].temp_factor = 0.0 ;

      /* Free some memory */
      for( j = 1 ; j <= Rotated_at_Origin_Mobile_Structure.length ; j ++ ) {
        free( Rotated_at_Origin_Mobile_Structure.Residue[j].Atom ) ;
        free( Translated_and_Rotated_at_Origin_Mobile_Structure.Residue[j].Atom ) ;
      }
      free( Rotated_at_Origin_Mobile_Structure.Residue ) ;
      free( Translated_and_Rotated_at_Origin_Mobile_Structure.Residue ) ;

    }

  }

  /* End of Loop */

  Centres_Structure.length = length ;

  /* Write pdb file */
  write_structure_to_pdb( Centres_Structure , "centres.pdb" ) ;

/************/

  /* Free the memory */

  for( i = 1 ; i <= Centres_Structure.length ; i ++ ) {
    free( Centres_Structure.Residue[i].Atom ) ;
  }
  free( Centres_Structure.Residue ) ;

  for( i = 1 ; i <= Origin_Mobile_Structure.length ; i ++ ) {
    free( Origin_Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Mobile_Structure.Residue ) ;

/************/

  printf( "Finished\n\n" ) ;

  return( 0 ) ;

} /* end main */
