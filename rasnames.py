#Python function to parse a mix of file names and dirctories
#HISTORY
#21Sep16 GIL initial version based on ras.py
import os
import glob

def splitNames( names, typea=".ast", typeb=".hot", doDebug=False):
    """
    splitNames parses the names in a list, looking for "typea" and "typeb".
    splitNames also expands directories
    to include all files in any directories found
    Inputs
    names - list of directories and files names
    typea - first type of file to search for, must not be blank
    typeb - second type of file to search for, may be blank (ignored)
    The function only searches one directory down from the current directory.
    """
    typea = str(typea)   # make sure we're matching strings
    typeb = str(typeb)
    if len(typea) < 1:
        print("splitNames: Invalid input type")
        return [], 0
    # flag whether searching for a second type of file
    findb = len(typeb) > 0
        
    if isinstance( names, str):
        names = names.split()
    nnames = len(names)
    count = 0
    tempnames = []
    #must join a wildcard match
    stara = "*"+typea
    starb = "*"+typeb
    if doDebug:
        print("Type A: '%s',  Type B: '%s'" % ( typea, typeb))

    for aname in names:
        
        if os.path.isdir(aname):  # if a directory name
            newnames = list(glob.glob(os.path.join(aname,stara)))
            if doDebug:
                print(newnames)
            if count == 0:
                tempnames = newnames
            else:
                for anothername in newnames:
                    tempnames.append(anothername)
            count = count + len(newnames)
            if findb:
                newnames = list(glob.glob(os.path.join(aname,starb)))
                if count == 0:
                    tempnames = newnames
                else:
                    for anothername in newnames:
                        tempnames.append(anothername)
                count = count + len(newnames)
        else:  # else this is just a file name
            if typea in aname:
                tempnames.append(aname)
                count = count + 1
            elif typeb in aname:
                tempnames.append(aname)
                count = count + 1
    # need files in time order which happens to be alphabetical order, too
    typenames = sorted(tempnames)
    if doDebug:
        print ("Found %d files name of types: '%s' '%s'" % (count, typea, typeb))
    # end of splitNames
    return typenames, count
        