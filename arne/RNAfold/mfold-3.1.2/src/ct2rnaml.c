const char COMMENT[]  = "Output from mFold";
const char DTD_DESC[] = "-//University of Montreal//RNAML 1.1//EN";
const char DTD_URL[]  = "http://www-lbit.iro.umontreal.ca/rnaml/current/rnaml.dtd";

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "getopt.h"

struct ctLine
{
  int seqnum;
  char base;
  int previous;
  int next;
  int paired;
  int historic;
  int strand;
};

struct ctHead
{
  int num_bases;
  int has_dG;
  int has_initial;
  float dG_val;
  float initial_val;
  char name[250];
  int is_rna, is_dna;
  int num_strands;
};

struct ssLine
{
  int seqnum;
  char base;
  double x;
  double y;
  int paired;
};

int readCtFile(FILE*, struct ctHead*, struct ctLine**);
int readSsFile(FILE*, struct ssLine*, int);
int isMatch(int, struct ctLine*, struct ssLine*);
int diff(int, struct ctLine*, struct ctLine*);
void makeValidId(char*);

const struct option OPTIONS[] = {
  {"melt", 1, 0, 'l'},
  {"method", 1, 0, 'm'},
  {0, 0, 0, 0}
};

int main(int argc, char** argv)
{
  char buffer[100];
  char *method, *melt;
  struct ctHead head, oldHead;
  int seqnew_bool, seq_num;
  int i, model_num, strand, strandBegin, strandEnd, helixLength;
  struct ctLine *lines, *oldLines;
  struct ssLine* ssLines;
  FILE *ctFile, *ssFile;
  FILE *output;
  
  if (argc == 1)
    {
      printf("Usage: %s [options] ct_file\n", argv[0]);
      printf("[-l, --melt=<melting temperature>]\n");
      printf("[-m, --method=<method text>]\n");
      return 1;
    }

  melt = NULL;
  method = NULL;
  while ((i = getopt_long(argc, argv, "l:m:", OPTIONS, 0)) != -1)
    if (i == 'l')
      melt = optarg;
    else if (i == 'm')
      method = optarg;

  if (optind >= argc)
    {
      fprintf(stderr, "Error: ct file not specified\n");
      return 1;
    }

  if (!strcmp(argv[optind] + strlen(argv[optind]) - 3, ".ct"))
    strcpy(argv[optind] + strlen(argv[optind]) - 3, "");

  strcpy(buffer, argv[optind]);
  strcat(buffer, ".ct");
  ctFile = fopen(buffer, "r");

  strcpy(buffer, argv[optind]);
  strcat(buffer, ".rnaml");
  output = fopen(buffer, "w");
  
  if (!ctFile || !output)
    {
      printf("An error occurred while opening the file\n");
      return 1;
    }
  
  fprintf(output, "<?xml version=\"1.0\"?>\n");
  fprintf(output, "<!DOCTYPE rnaml PUBLIC \"%s\" \"%s\">\n\n", DTD_DESC, DTD_URL);
  fprintf(output, "<rnaml version=\"1.1\" comment=\"%s\">\n\n", COMMENT);
  fprintf(output, "  <analysis id=\"mfold\">\n");
  fprintf(output, "    <program>\n");
  fprintf(output, "      <prog-name>mFold</prog-name>\n");
  fprintf(output, "      <prog-version>3.1</prog-version>\n");
  fprintf(output, "    </program>\n");
  fprintf(output, "  </analysis>\n\n");
  
  seqnew_bool = 1;
  seq_num = 1;
  model_num = 0;
  lines = NULL;
  ssLines = NULL;
  
  while (1)
    {
      bcopy(&head, &oldHead, sizeof(struct ctHead));
      bzero(&head, sizeof(struct ctHead));
      oldLines = lines;

      if (!readCtFile(ctFile, &head, &lines))
	break;

      makeValidId(head.name);

      if ((seq_num != 1 || model_num != 0) && (head.num_bases != oldHead.num_bases || strcmp(head.name, oldHead.name) || diff(head.num_bases, lines, oldLines)))
	{
	  seqnew_bool = 1;
	  ++seq_num;
	  model_num = 1;
	}
      else
	++model_num;

      strcpy(buffer, argv[optind]);
      strcat(buffer, ".ss");
      ssFile = fopen(buffer, "r");
      if (!ssFile)
	{	
	  sprintf(buffer, "%s_%d.ss", argv[optind], model_num);
	  ssFile = fopen(buffer, "r");
	}
      if (ssFile)
	{
	  free(ssLines);
	  ssLines = (struct ssLine*) malloc(head.num_bases * sizeof(struct ssLine));
	  if (!readSsFile(ssFile, ssLines, head.num_bases))
	    {
	      free(ssLines);
	      ssLines = NULL;
	    }
	  if (!isMatch(head.num_bases, lines, ssLines))
	    {
	      fprintf(stderr, "Warning: %s is inconsistent and will be ignored\n", buffer);
	      free(ssLines);
	      ssLines = NULL;
	    }
	}
      
      
      /*this is the beginning of the output to the file*/
      if (seqnew_bool)
	{  
	  if (seq_num != 1)
	    {
	      fprintf(output, "    </structure>\n");
	      fprintf(output, "  </molecule>\n");
	    }
	  
	  if ((seq_num % 10) == 1)
	    fprintf(output, "\n  <!-- %dst sequence and its associated foldings ******************* -->\n", seq_num);
	  else if ((seq_num % 10) == 2)
	    fprintf(output, "\n  <!-- %dnd sequence and its associated foldings ******************* -->\n", seq_num);
	  else if ((seq_num % 10) == 3)
	    fprintf(output, "\n  <!-- %drd sequence and its associated foldings ******************* -->\n", seq_num);
	  else
	    fprintf(output, "\n  <!-- %dth sequence and its associated foldings ******************* -->\n", seq_num);
	  
	  if (head.is_dna && head.is_rna)
	    fprintf(output, "  <!-- Warning: type unknown as sequence contains both Ts and Us -->\n");
	  if (head.is_dna && !head.is_rna)
	    fprintf(output, "  <molecule id=\"%s\" type=\"dna\">\n", head.name);
	  else if (head.is_rna && !head.is_dna)
	    fprintf(output, "  <molecule id=\"%s\" type=\"rna\">\n", head.name);
	  else
	    fprintf(output, "  <molecule id=\"%s\">\n", head.name);

	  for (strand = 1; strand <= head.num_strands; ++strand)
	    {
	      for (strandBegin = 0; lines[strandBegin].strand < strand; ++strandBegin);
	      for (strandEnd = strandBegin; strandEnd < head.num_bases && lines[strandEnd].strand == strand; ++strandEnd);

	      if (lines[strandBegin].previous == lines[strandEnd - 1].seqnum && 
		  lines[strandEnd - 1].next == lines[strandBegin].seqnum)
		fprintf(output, "    <sequence circular=\"true\" strand=\"%d\">\n", strand);
	      else
		fprintf(output, "    <sequence circular=\"false\" strand=\"%d\">\n", strand);
	      fprintf(output, "      <seq-data>\n        ");

	      for (i = 0; strandBegin + i < head.num_bases && i < strandEnd; ++i)
		if (lines[strandBegin + i].strand == strand)
		  {
		    fprintf(output, "%c", lines[strandBegin + i].base);
		    if (i % 60 == 59)
		      fprintf(output, "\n        ");
		    else if (i % 10 == 9)
		      fprintf(output, " ");
		  }
	      fprintf(output, "\n      </seq-data>\n");
	      fprintf(output, "    </sequence>\n\n");
	    }
	  fprintf(output, "    <structure analysis-ids=\"mfold\">\n");
	  seqnew_bool = 0;
	}

      fprintf(output, "      <!-- Start of folding %d for sequence %d ************************ -->\n", model_num, seq_num);
      if (melt)
	fprintf(output, "      <model id=\"_%d.%d\" comment=\"Melting temperature: %s\">\n", seq_num, model_num, melt);
      else
	fprintf(output, "      <model id=\"_%d.%d\">\n", seq_num, model_num);
      fprintf(output, "        <model-info>\n");
      if (method)
	fprintf(output, "          <method>%s</method>\n", method);
      
      if (head.has_dG)
	{
	  if (head.has_initial)
	    fprintf(output, "          <free-energy comment=\"refined\">%.2f</free-energy>\n", head.dG_val);
	  else
	    fprintf(output, "          <free-energy>%.2f</free-energy>\n", head.dG_val);
	}      
      if (head.has_initial)
	fprintf(output, "          <free-energy comment=\"initial\">%.2f</free-energy>\n", head.initial_val);
      
      fprintf(output, "        </model-info>\n\n");
      
      /*if (model_num == 1)
	{
	  for (i = 0; i < head.num_bases; i++)
	    {
	      fprintf(output, "        <base>\n");
	      if (head.num_strands > 1)
		fprintf(output, "          <strand>%d</strand>\n", lines[i].strand);
	      fprintf(output, "          <position>%d</position>\n", lines[i].historic);
	      fprintf(output, "          <base-type>%c</base-type>\n", lines[i].base);
	      fprintf(output, "        </base>\n");
	    }
	  fprintf(output, "\n");
	  }*/
      
      fprintf(output, "        <str-annotation>\n");
      for (i = 0; i < head.num_bases; i++)
	if (lines[i].seqnum < lines[i].paired)
	  {
	    helixLength = 1;
	    if (lines[i].paired > 1)
	      while (i + helixLength < head.num_bases && lines[i + helixLength].paired == lines[i].paired - helixLength)
		++helixLength;

	    if (helixLength > 1)
	      fprintf(output, "          <helix>\n");
	    else
	      fprintf(output, "          <base-pair>\n");
	    fprintf(output, "            <base-id-5p>\n");
	    fprintf(output, "              <base-id>\n");
	    if (head.num_strands > 1)
	      fprintf(output, "                <strand>%d</strand>\n", lines[i].strand);
	    fprintf(output, "                <position>%d</position>\n", lines[i].historic) ;
	    fprintf(output, "              </base-id>\n");
	    fprintf(output, "            </base-id-5p>\n");
	    fprintf(output, "            <base-id-3p>\n");
	    fprintf(output, "              <base-id>\n");
	    if (head.num_strands > 1)
	      fprintf(output, "                <strand>%d</strand>\n", lines[lines[i].paired - 1].strand);
	    fprintf(output, "                <position>%d</position>\n", lines[lines[i].paired - 1].historic);
	    fprintf(output, "              </base-id>\n");
	    fprintf(output, "            </base-id-3p>\n");
	    if (helixLength > 1)
	      {
		fprintf(output, "            <length>%d</length>\n", helixLength);
		fprintf(output, "          </helix>\n");  
		i += (helixLength - 1);
	      }
	    else
	      fprintf(output, "          </base-pair>\n");  
	  }
      fprintf(output, "        </str-annotation>\n");
      if (ssLines)
	{
	  fprintf(output, "\n        <secondary-structure-display>\n");
	  for (i = 0; i < head.num_bases; ++i)
	    {
	      fprintf(output, "          <ss-base-coord>\n");
	      fprintf(output, "            <base-id>\n");
	      if (head.num_strands > 1)
		fprintf(output, "              <strand>%d</strand>\n", lines[i].strand);
	      fprintf(output, "              <position>%d</position>\n", lines[i].historic) ;
	      fprintf(output, "            </base-id>\n");
	      fprintf(output, "            <coordinates>%lg %lg</coordinates>\n", ssLines[i].x, ssLines[i].y);
	      fprintf(output, "          </ss-base-coord>\n");
	    }
	  fprintf(output, "        </secondary-structure-display>\n");
	}
      fprintf(output, "      </model>\n");
      fprintf(output, "      <!-- End of folding %d for sequence %d ************************ -->\n\n", model_num, seq_num);
      free(oldLines);
    }	
  
  fprintf(output, "    </structure>\n");
  fprintf(output, "  </molecule>\n\n");
  fprintf(output, "</rnaml>\n");
  
  return 0;
}

int readCtFile(FILE* file, struct ctHead* head, struct ctLine** lines)
{
  int i;
  int currentStrand = 1;
  char current_line[300];

  if (!fgets(current_line, 300, file))
    return 0;

  if (sscanf(current_line, "%d dG = %G [initially %g] %[^\n]",
	     &head->num_bases, &head->dG_val, &head->initial_val, head->name) == 4)
    {
      head->has_dG = 1;
      head->has_initial = 1;
    }
  else if (sscanf(current_line, "%d dG = %g %[^\n]", &head->num_bases, &head->dG_val, head->name) == 3)
    head->has_dG = 1;
  else if (sscanf(current_line, "%d ENERGY = %g [initially %g] %[^\n]",
		  &head->num_bases, &head->dG_val, &head->initial_val, head->name) == 4)
    {
      head->has_dG = 2;
      head->has_initial = 1;
    }
  else if (sscanf(current_line, "%d ENERGY = %g %[^\n]", &head->num_bases, &head->dG_val, head->name) == 3)
    head->has_dG = 2;
  else if (sscanf(current_line, "%d %[^\n]", &head->num_bases, head->name) == 2)
    ;
  else
    return 0;

  if (head->num_bases <= 0)
    return 0;

  while (head->name[strlen(head->name) - 1] == ' ')
    head->name[strlen(head->name) - 1] = '\0';

  if (!(*lines = (struct ctLine*) malloc(sizeof(struct ctLine) * head->num_bases)))
    return 0;

  for (i = 0; i < head->num_bases; i++)
    {
      if (!fgets(current_line, 100, file))
	return 0;
      if (sscanf(current_line, "%d %c %d %d %d %d", &(*lines)[i].seqnum, &(*lines)[i].base, &(*lines)[i].previous, &(*lines)[i].next, &(*lines)[i].paired, &(*lines)[i].historic) != 6)
	return 0;

      if ((*lines)[i].seqnum != i + 1)
	{
	  printf("The sequence does not go up in sequentially increasing numeric order\n");
	  return 0;
	}

      if (i > 0 && (*lines)[i].previous != i && (*lines)[i - 1].next != i + 1)
	++currentStrand;

      (*lines)[i].strand = currentStrand;

      if (toupper((*lines)[i].base) == 'T')
	head->is_dna = 1;
      else if (toupper((*lines)[i].base) == 'U')
	head->is_rna = 1;
    }

  head->num_strands = currentStrand;

  return 1;
}

int readSsFile(FILE* file, struct ssLine* lines, int n)
{
  int i;
  char current_line[300];

  for (i = 0; i < n; ++i)
    {
      if (!fgets(current_line, 300, file))
	return 0;
      if (sscanf(current_line, "%d %c %lg %lg %*d %d", &lines[i].seqnum, &lines[i].base, &lines[i].x, &lines[i].y, &lines[i].paired) != 5)
	return 0;
      if (lines[i].seqnum != i + 1)
	return 0;
    }

  return 1;
}

int isMatch(int n, struct ctLine* ctLines, struct ssLine* ssLines)
{
  int i;

  for (i = 0; i < n; ++i)
    if (ctLines[i].base != ssLines[i].base || ctLines[i].paired != ssLines[i].paired)
      return 0;

  return 1;
}

int diff(int size, struct ctLine* new, struct ctLine* old)
{
  int i;

  for (i = 0; i < size; ++i)
    if (new[i].base != old[i].base)
      return 1;

  return 0;
}

void makeValidId(char* str)
{
  char buffer[250];

  if (('0' <= *str && *str <= '9') || *str == '-' || *str == '.')
    {
      strcpy(buffer, str);
      strcpy(str, "_");
      strcat(str, buffer);
    }

  for (; *str; ++str)
    if (!isalnum(*str) && *str != '_' && *str != '-' && *str != '.')
      *str = '_';
}
