## Relational Databases >>

For one-off projects it is probably okay to work directly with raw data from spreadsheets or delimited files, because the act of creating anything more sophisticated to hold the raw data will take longer than just performing the analysis straight-away. For long-term projects, spanning years, or sometimes months, I STRONGLY recommend creating some sort of database to receive and standardize raw data inputs before beginning analysis. In these cases, investing a day or two to set-up the database can save an enormous amount of time down the road. Also, it will make it easier for someone to pick-up where others have left off. I've used two database systems worth mentioning, though there are many to choose from, and I'll briefly mention their pros and cons here...

### Microsoft Sharepoint online:

If you are a PC user and appreciate Microsoft products, then you may want to take a look at Sharepoint Online. Sharepoint can do many things which I have not fully explored, however, I am familiar with its relational database structure. If you don't know what a relational database is, Google it. Anyone with a microsoft account can purchase an affordable subscription to sharepoint (along the order of 5 USD per month). Once you have a Sharepoint site, you can make it public or private, and you can start to set-up linked data tables to hold all of your datasets. 

#### Here are some nice features of a relational database via Sharepoint online:
* You can easily give access to anyone else that has a microsoft account, and you can control which objects within the database they can view, edit, create, delete, etc.
* Sharepoint syncs easily with the Microsoft Access desktop app. If you know how to use Access, well now you can work with Access but your data is not stored locally, it is stored on the Sharepoint server. As such, anyone with Access can work with your database from anywhere in the world
* Syncing with Access has some really big advantages, namely that you can create data entry forms for people to work with. Access has some nifty ways of integrating subforms and lookup columns, which are not nearly as create on your own.
* You can create public views and data entry forms online, which can be really useful

#### Here are the bad features of a relational database via Sharepoint online:
* Above all other things, I find that Microsoft Access is glitchy, especially in so far as granting access to an online database. There are security requirements that need to be met that are not easy to figure out when you encounter an error
* Access only runs on Windows. For me, I have had to create virtual machines to use Access, which can be a nuisance. The more steps there are to run an application the more likely it is that they are interrupted when Microsoft, Apple, or some other third party updates its software.

### SQL:

Create a database with structured query language (SQL). There are a lot of software products that will help you do this, e.g. [mysql](https://www.mysql.com/) or [postgresql](https://www.postgresql.org/). I currently use postgresql. Arguably, the learning curve is a little steeper here, but once you understand how to set-up and work with a relational database through sql the sky is the limit. 
#### Here are some obvious advantages:
* There are endless online resources for learning to use SQL
* The language is intuitive, and actually it is what Sharepoint and Access are using behind the scenes, so if you can execute the queries without the help of a GUI, why not.
* Easy to give others access to your database, but perhaps harder for them to work with it, unless they also know SQL
* Can be stored and accessed remotely. For example, I put my postgresql database onto [Amazon Web Services](https://aws.amazon.com/free), which has some vary affordable hosting services. AWS also has most relation database platforms ready to launch, so you just need to open an account, restore an SQL dump file, and you are ready to work.

#### Here are some disadvantages:
* If working from the command line or raw code makes you nervous, you would have to get over that
* There is no one simple solution for making your database accessible to the lay person. But with determination anything will be possible
* Anytime you do things on the command line it is a little easier to make catastrophic mistakes. So you need to be familiar with how to create snapshots or backups of your database.

Whatever you choose to do in the end, you will save an enormous amount of time inputing raw data into relational databases first, and then extracting what you need for a particular analysis. Relational databases force your data to have referential integrity, so it is a little bit harder for someone to screw up an entire analysis with a simple typo.



## Questions or Comments
Please get in touch through my webpage [GideonErkenswick.com](https://gideonerkenswick.com/contact/)
