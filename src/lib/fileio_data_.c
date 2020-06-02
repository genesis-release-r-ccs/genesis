
#include <stdio.h>

#define MaxFile  10

static FILE* Files [MaxFile] = { 0 };

void fd_open__(char* filename, const int* read_only, int* unit_no,
			  int fnlen)
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
		Files[i] = fopen(filename, (*read_only) == 1 ? "r" : "w");
		break;
	  }

  if (i == MaxFile)
	*unit_no = -1;
}

void fd_close__(int* unit_no)
{
  if (*unit_no < 0 || *unit_no >= MaxFile)
	return ;
  fclose(Files[*unit_no]);
  Files[*unit_no] = 0;
}

void fd_flen__(int* unit_no, long long* flen)
{
  FILE* fp = Files[*unit_no];
  fseeko64(fp, 0, SEEK_END);
  *flen = ftello64(fp);
  fseeko64(fp, 0, SEEK_SET);
}

void fd_read__(int* unit_no, long long *pos, unsigned char* b, long long *blen)
{
  FILE* fp = Files[*unit_no];
  fseeko64(fp, (*pos)-1, SEEK_SET);
  fread(b, *blen, 1, fp);
}

void fd_write__(int* unit_no, unsigned char* b, long long *blen)
{
  FILE* fp = Files[*unit_no];
  fwrite(b, *blen, 1, fp);
}
