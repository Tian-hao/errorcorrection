#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
	double counta,countc,countt,countg,deletion,insertion;
	double count_total;

	FILE  *nuc ;
	FILE  *fre ;

	char ch1,ch3;
	char ch2[999];
	double m;
	double x;
	double transition,transversion;
	
	char *insert;
	int ins_pos;

	FILE *output;
        int i;
	
	char *inputname;
	char *outputname;
	char *tempname;

        inputname = (char*)malloc(sizeof(char)*50);
        outputname = (char*)malloc(sizeof(char)*50);
	tempname = (char*)malloc(sizeof(char)*50);
	//scanf("%s",inputname);
	sscanf(argv[2],"%s",inputname);
        printf("good memory\n");

        sscanf(inputname,"../alignment/%s",tempname);
        sprintf(outputname,"../graph/%s",tempname);

        nuc = fopen(inputname,"r");
        fre = fopen(inputname,"r");
        output = fopen(outputname,"w");

        printf("reading %s\n",inputname);
        printf("writing %s\n",outputname);
        for(i=0;i<88;i++)
        {
            fseek(nuc,i,SEEK_SET);
            fseek(fre,89,SEEK_SET);
            counta = 0; countt = 0; countc = 0; countg = 0; deletion = 0; insertion = 0;

            while(1)
            {
                ch1 = fgetc(nuc);
                if (ch1 == EOF) break;
                fscanf(fre,"\t%lf\t%lf",&m,&x);
		insert = (char*)malloc(sizeof(char)*50);
		fscanf(fre,"\t%s",insert);
		ins_pos = -1;
		while(*insert)
		{
			sscanf(insert,"%d",&ins_pos);
			if(ins_pos == i+1)
			{
				insertion += m;
				break;
			}
			else
			{
				insert += 4;
				continue;
			}
		}

                if (ch1 == 'A') counta += m;
                else if (ch1 == 'C') countc += m;
                else if (ch1 == 'T') countt += m;
                else if (ch1 == 'G') countg += m;
		else if (ch1 == 'N') 
		{
			counta += 0.25*m;
			countc += 0.25*m;
			countt += 0.25*m;
			countg += 0.25*m;
		}
		else if (ch1 == '-') deletion += m;
                fgets(ch2,999,nuc);
                fgets(ch2,999,fre);
                fseek(nuc,i,SEEK_CUR);
                fseek(fre,89,SEEK_CUR);
            }
            fprintf(output,"%d\t",i+1);
            count_total = counta + countt + countc + countg + deletion + insertion;
	    counta = counta/count_total;
	    countt = countt/count_total;
	    countc = countc/count_total;
	    countg = countg/count_total;
	    if(counta>0.5){ transition = countg; transversion = countc + countt;}
	    if(countt>0.5){ transition = countc; transversion = counta + countg;}
	    if(countc>0.5){ transition = countt; transversion = counta + countg;}
	    if(countg>0.5){ transition = counta; transversion = countc + countt;}
            fprintf(output,"%.10lf\t",transition);
            fprintf(output,"%.10lf\t",transversion);
	    x = insertion/count_total;
	    fprintf(output,"%.10lf\t",x);
	    x = deletion/count_total;
	    fprintf(output,"%.10lf\n",x);
        }

        fclose(nuc);
        fclose(fre);
        fclose(output);
	
	return 0;
}
