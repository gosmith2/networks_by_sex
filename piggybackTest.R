##Piggyback fiddling

library(piggyback)
Sys.setenv(GITHUB_TOKEN="7423a37de1b5d3913210273a64d0d89e65e5c8f9")

#For more info throughout this, see
#https://cran.r-project.org/web/packages/piggyback/vignettes/intro.html

##-----First upload-----##

#If there are no releases on the repository, pb_new_release()
#will return an error. Create a release on github itself.
#For this, I named it "0.01"

path<-"/Users/thebo/Documents/UCR_postdoc/networks_by_sex/piggyTest.csv"
  #Or whatever path you need for the file. should work with
  #a url?

pb_upload(path,
          tag="0.01") 
#Tag is the name of the release you created on Github.
#pb_upload optionally takes a repository path if you aren't
#already in the correct repository

#the file is now associated with the repository, and should
#stay as is. 

#so to test things, I deleted the original file from my computer

##-----Using the files in the release------#

#You'll need a token for this... I don't know as much about
#them, but I made one and set it temporarily as an environment
#variable that I then took out of the script

pb_download("piggyTest.csv",
            tag="0.01")

#You can set things like the destination if you like; it 
#defaults to your local repo folder. 

PT<-read.csv("piggyTest.csv")
PT


#At this point I changed the original csv file on my computer
#to test things, adding a 6th line of "data"

pb_download("piggyTest.csv",
            tag="0.01")
#this time dowload didn't seem to do much: just said that 
#the files are up to date. The file on my computer still 
#has that 6th line. So this probably only overwrites your
#stuff if its more recent, which should almost never happen?


#As a further test, I created a new release and uploaded
#the new 6-line test csv to it. 
pb_new_release(tag="v0.01")

pb_upload("piggyTest.csv",tag="v0.01")

#at this point I can get either the 5 or 6-line versions
#depending on which release I download from. 

#Last thing: I added a 7th line, and tried to upload it
#to the new release

pb_upload("piggyTest.csv",
          tag="v0.01")

#deleted the version on my computer, downloaded the most
#recent release:

pb_download()

PT1<-read.csv("piggyTest.csv")
PT1

#has all 7 rows. So you "push" by just re-uploading
#something with the same file name to the release
