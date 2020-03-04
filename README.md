#shinyNAP
#Setting up ShinyNAP
#### Using Ubuntu 18.04 on Amazon AWS EC2
The purpose of this markdown file is to reproduce the R (Shiny) application **shinyNAP**, or shiny NeoAntigenPortal deployed onto an Amazon AWS EC2 server instance (or other LINUX/UNIX-based systems). ShinyNAP is an application designed to provide oncology researchers with a visualization tool for protein immunogenicity. 

Using AWS EC2
---------


For general setup using Amazon AWS: 

- Log into [Amazon] (https://signin.aws.amazon.com/signin?redirect_uri=https%3A%2F%2Fus-east-2.console.aws.amazon.com%2Fec2%2Fv2%2Fhome%3Fregion%3Dus-east-2%26state%3DhashArgs%2523Instances%253Asort%253DkeyName%26isauthcode%3Dtrue&client_id=arn%3Aaws%3Aiam%3A%3A015428540659%3Auser%2Fec2&forceMobileApp=0)
    1. Click **create instance**
    2. Click **launch screen**
    3. Select **Ubuntu18.04** (free instance)
    4. Select **General Purpose - t2.micro**
    5. Click **Review and Launch** in the bottom-left

- Change security settings for instance
    1. in the SSH row, change source to **My IP**
    2. Click **Add Rule** to add custom TCP rule
    3. Under **Port Range** enter **3838** (port for R shiny server)
    4. Under **Port Range** enter **8787** (port for Rstudio Server)
    5. Change **Source** to **Anywhere**
    6. Click **Review and Launch** in the bottom left
    7. Click **Launch** to get the private key to SSH into the EC2 instance
    8. Create new key pair **NAPkey** to download **NAPkey.pem** and save to an easily accessible directory
        
        - Secure the key by using ```chmod 400 NAPkey.pem```
    9. Copy address to instance under **Public DNS (IPv4)**
    
        -Ex. *ec2-18-223-112-142.us-east-2.compute.amazonaws.com*
        
####Connect using SSH
In AWS console, select ***connect***, and copy the path in the dialogue box with the ***SAME*** directory the Key file is in:
```sudo ssh -i "NAPkey.pem" ubuntu@ec2-18-223-112-142.us-east-2.compute.amazonaws.com```

###Install R, etc. on server
-----
You have now successfully logged into the server instance you created!!! Let's use it. The following instructions/commands are adapted from [Kimberly Coffee] (http://www.kimberlycoffey.com/blog/2016/2/13/mlz90wjw0k76446xkg262prvjp0l8u). You can also follow a useful video from [Youtube] (https://www.youtube.com/watch?v=hglgkFfRqyQ) that also uses Kimberly's commands.

```sudo su -c "echo 'deb http://archive.linux.duke.edu/cran/bin/linux/ubuntu trusty/' >> /etc/apt/sources.list"```

```sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9```

```sudo apt-get update```

```sudo apt-get upgrade -y #when receiving warning here, click continue!```

```udo apt-get dist-upgrade -y ```

```sudo apt-get install r-base -y ```

Install R packages here:

```sudo su -c "R -e \"install.packages('shiny', repos='https://cran.rstudio.com/')\"" ```

```sudo su -c "R -e \"install.packages('rmarkdown', repos='https://cran.rstudio.com/')\"" ```

``` sudo su -c "R -e \"install.packages('ggplot2', repos='https://cran.rstudio.com/')\""```

``` sudo su -c "R -e \"install.packages('tidyverse', repos='https://cran.rstudio.com/')\""```

``` sudo su -c "R -e \"install.packages('plotly', repos='https://cran.rstudio.com/')\""```

```sudo su -c "R -e \"install.packages('doMC', repos='https://cran.rstudio.com/')\"" ```

For other packages...
```sudo su -c "R -e \"install.packages('PACKAGE NAME HERE', repos='https://cran.rstudio.com/')\"" ```

####Install shiny server
Use the most-updated shiny server [link] (https://www.rstudio.com/products/shiny/download-server/) and use the following commands in the server.
```sudo apt-get install gdebi-core```

```sudo gdebi shiny-server-1.5.9.923-amd64.deb```

Test the above by using your server's IP address (found in the console)

<http://18.223.112.142/3838>

####Install Rstudio server
Use the most-updated Rstudio server [link] (https://www.rstudio.com/products/rstudio/download-server/)

```wget https://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb```

```sudo gdebi rstudio-server-1.1.463-amd64.deb```

####Add user for login to Rstudio
You can obviously choose your own login. It will also prompt you to enter your own password...
```sudo adduser shinynapadmin```

####Add sudo user priveleges
```sudo vim /etc/sudoers```

Scroll down to **allow members of group sudo to execute any command. Press "i" to insert:

```%shinynapadmin ALL=(ALL:ALL) ALL```

Check sudoers:

```sudo cat /etc/sudoers```

##Transfer NetMHC files from LOCAL to server
In the folder where the key is located, enter the following command to transer R scripts to the server (you can do this, or git clone from the repo. ShinyNAP should already have a configured version of NetMHC included in the app.

```scp -i NAPkey.pem *.R```

Now, download netMHC to a local location. You will be transferring the unlocked version from the DTU Bioinformatics [website] (http://www.cbs.dtu.dk/services/NetMHC/). You will be prompted to enter a **.edu** email address in order to download.

Once you have the netMHC files, transfer them to the server using the ubuntu/server address from your EC2 console. ```sudo scp -i NAPkey.pem -r netMHC-4.0  ubuntu@ec2-3-16-109-251.us-east-2.compute.amazonaws.com:~```


####Configure NetMHC software on server

Install tsch to execute netMHC:

```sudo apt-get install csh```

```sudo apt-get install tcsh```

Set shortcut to "data" folder for netMHC. In the netMHC-4.0 directory, change the path from its default to **/home/ubuntu/netmhc_files/netMHC-4.0** using the following command

```sudo pico netMHC```

Download "test" data from server using the link provided in netMHC's download [email] (http://www.cbs.dtu.dk/services/NetMHC-4.0/data.tar.gz).

Transfer gzipped files to server:

```sudo scp -i NAPkey.pem -r data.tar.gz  ubuntu@ec2-3-16-109-251.us-east-2.compute.amazonaws.com:~/netMHC-4.0/```

Unzip the files:

```gunzip -c data | tar xvf - ```

Test netMHC per instructions:

```sudo ../netMHC test.fsa > test.fsa.myout```

```sudo ../netMHC -p test.pep > test.pep.myout```

Once netMHC is configured, move (or copy) netMHC to the www directory in the shinyNAP app directory www (netMHC should already be configured when the repo is cloned).

####Restart server after installation of all updates

```sudo reboot```

```sudo apt-get update```

#Clone ShinyNAP git repo
####in Rstudio
Before we get started, enable sudoer permissions for your usernames/shiny's usernames so shiny can read/write/execute.
```sudo pico /etc/sudoers```

- add **%shinynapadmin ALL=(ALL:ALL) ALL**
- add **%ubuntu ALL=(ALL:ALL) ALL**
- add **%shiny ALL=(ALL:ALL) ALL**

The moment you've all been waiting for...  

Once you're logged into Rstudio by entering your.IP.address/8787:

1. Click **new project**
2. Click **version control**
3. Click **git**
4. Enter address **https://github.com/miosisoniii/shinyNAP**
5. Enter your git username
6. Enter password

The app ***should*** run when using Rstudio on the server to run the app, ***with*** NetMHC.



If there are any questions please email at miosisoniii@gmail.com or slack me at miosisoniii.