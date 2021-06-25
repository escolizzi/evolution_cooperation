# evolution_cooperation
Source code for the results in the article: Colizzi, Hogeweg 2016. High cost enhances cooperation through the interplay between evolution and self-organisation

##############################################################################################################

  In this directory you can find the code we used to obtain results reported in the publication:
  "High cost enhances cooperation through the interplay between evolution and self-organisation",
  E.S. Colizzi & P. Hogeweg, 2016. BMC Evolutionary Biology.
  
  The code is written in c, and makes use of the CASH libraries, 
  which implement Cellular Automata dynamics, and are also freely available at the link:
  http://bioinformatics.bio.uu.nl/rdb/software.html
  
  For the sake of simplicity, I made large use of the environment provided by the student-oriented demo of CASH
  (which is just a thin wrapper around the libraries that makes them immediately usable),
  which is available as a tar ball (cash2-s.[year].tar.gz) at the link:
  http://bioinformatics.bio.uu.nl/BINF/Programs_used_in_course/
  
  In the same page there is a manual on how to use CASH, and my code is usable just like any other file.
  Libraries and code are compiled with gcc under Ubuntu, and I guess that most of Linux systems should work with hardly any hacking.
  I do not know whether anything here compiles under any other operating system.
  
  You can freely download and modify my code.
  My code uses the Unlicense, so do with it what you want!
  Citing us if you publish anything with it 'is' highly appreciated, though :P

  Once you compiled the code, type something like: "./public_good_5 jksdfnksjd" to see the command line options.
  The program generates three types of outputs:
  - "dataPB.txt" is a text file which contains random samples of the entire population at regular time intervals.
    Beware: This file may get big if you sample too many individuals too frequently!
  - "moviePB" is a directory which contains a collection of numbered png files that can be wathed subsequently very fast to produce the illusion of a movie.
  - "backupPB" is a collection of text files each containing a whole dump of the information of the lattice in ordered fashion. 
  
  Do tell me about any bug you may find.
  ENJOY!
  Enrico Sandro Colizzi (email: istidina -at- gmail -dot- com)
  
##############################################################################################################
