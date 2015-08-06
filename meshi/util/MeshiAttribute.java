package meshi.util;
/*
 * All the additional attributes of Distance are type of DistanceAtribute
 */
public interface MeshiAttribute{
    final static int LENNARD_JONES_ELEMENT_ATTRIBUTE = 0;
    final static int EXCLUDED_VOLUME_ELEMENT_ATTRIBUTE = 1;
    final static int HYDROGEN_BONDS_ATTRIBUTE = 2;
    final static int CN_ATTRIBUTE = 3;
    final static int SEQUENCE_ALIGNMENT_COLUMN_ATTRIBUTE = 4;
    final static int SOLVATE_ALL_ATOM_ATTRIBUTE = 5;
    final static int SOLVATE_CA_ATTRIBUTE = 6;
    final static int SOLVATE_ROT1_ATTRIBUTE = 7;
    final static int EV_ROT1_ATTRIBUTE = 8;
    final static int SOLVATE_EXTRACTION_ATTRIBUTE = 9;    
    int key();
}
