/******************************************************************/
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <unistd.h>

#define MaxFile  128

static FILE* Files [MaxFile] = { 0 };

/******************************************************************/

void pr_fill(FILE* fp, long long *blen)
{
#define Unit 1024LL

  long long i;

  const char Data[Unit] = { };

  long long num = *blen / Unit;
  long long tail = *blen % Unit;

  for(i = 0; i < num; ++i)
	fwrite(Data, Unit, 1, fp);
  for(i = 0; i < tail; ++i)
	fwrite(&Data[0], 1, 1, fp);
}

/******************************************************************/

void pr_open_(char* filename, long long *blen, int* unit_no,
			  int fnlen, int mlen)
{
  int i;
  for(i = fnlen-1; i >= 0; --i)
	if (filename[i] == ' ')
	  filename[i] = 0;
	else
	  break;
	  
  for(i = 0; i < MaxFile; ++i)
	if (Files[i] == 0)
	  {
		*unit_no = i;

		char str[256];
		sprintf(str, "touch %s", filename);
		system(str);

		Files[i] = fopen(filename, "rb+");
		if (Files[i] == 0)
		  {
			printf("Pr_Open_> File open error [%s]\n", filename);
			*unit_no = -1;
			return ;
		  }

		pr_fill(Files[i], blen);

		//printf("OPENED : unit_no:%d  %s\n", *unit_no, filename);
		break;
	  }

  if (i == MaxFile)
	*unit_no = -1;
}

/******************************************************************/

void pr_close_(int* unit_no)
{
  if (*unit_no < 0 || *unit_no >= MaxFile)
	return ;
  fclose(Files[*unit_no]);
  Files[*unit_no] = 0;
}

/******************************************************************/

void pr_flen_(int* unit_no, long long* flen)
{
  FILE* fp = Files[*unit_no];
  fseeko(fp, (off_t)0, SEEK_END);
  *flen = (long long)ftello(fp);
}

/******************************************************************/

void pr_read_(int* unit_no, long long *pos, unsigned char* b, int *blen)
{
  FILE* fp = Files[*unit_no];
  fseeko(fp, (off_t)*pos, SEEK_SET);
  fread(b, *blen, 1, fp);
}

/******************************************************************/

void pr_write_(int* unit_no, long long *pos, unsigned char* b, int *blen)
{
  FILE* fp = Files[*unit_no];
  fseeko(fp, (off_t)*pos, SEEK_SET);
  fwrite(b, *blen, 1, fp);
}

/******************************************************************/
