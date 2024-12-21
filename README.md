## PCRpanel 
## Custom Amplicon Panels Designer is a tool for Amplicon Sequencing Technology, an ultra-high multiplex tiling PCR sequencing approach targeted DNA sequencing technology for next-generation sequencing (NGS) and nanopore sequencing (ONT).

By Ruslan Kalendar 

email: ruslan.kalendar@helsinki.fi

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 23 or higher

[Java Downloads](https://www.oracle.com/java/technologies/downloads/)


How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)



To run the project from the command line. Command-line options, separated by spaces. 
The executive file ```PCRpanel.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```
java -jar <PCRpanelPath>\dist\PCRpanel.jar <PCRpanelPath>\test\config.file

java -jar C:\PCRpanel\dist\PCRpanel.jar C:\PCRpanel\test\config.file 
```

### Basic usage: 
To enter parameters and specify the location of the target files and primer's file, you must specify this via a file on the command line. An example of such a file here (file name or extension does not matter):

> **config.file**
```
target_path=C:\PCRpanel\test\NG_008690.txt
target_path=C:\PCRpanel\test\NG_011731.txt
target_path=C:\PCRpanel\test\NG_013019.txt
target_path=C:\PCRpanel\test\NG_008847.txt
target_path=C:\PCRpanel\test\NC_000002.txt

target_primers=C:\PCRpanel\test\primers.txt

minPCR=50
maxPCR=500
minLen=18
maxLen=24
minTm=60
maxTm=62
minLC=78
3end=w
5end=
forwardtail=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
reversetail=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

```
## Specifies the coordinates of exons to be analyzed:
In the Genbank file, you must replace "mRNA" with "Panel" to indicate exon coordinates or any other target fragments to the program.

```
Panel            join(1..519,5421..5620,11239..11503,13230..13465,
                     34397..34557,40011..40143,43883..44077,45866..45984,
                     57668..58629)
```
