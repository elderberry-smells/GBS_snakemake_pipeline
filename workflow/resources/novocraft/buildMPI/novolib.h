/*
 * File:   novolib.h
 * Author: Colin Hercus
 * Copyright: 2012 Novocraft Technologies Sdn Bhd
 *
 * Created on March 1st, 2012, 11:57 AM
 */

#ifndef NOVOLIB_H
#define	NOVOLIB_H

#include <string>

    /*
     * Redirect stderr to /dev/null if it is same file as stdout
     */
    void chkstderr();
    /*
     * Check folder permisisons
     */
    bool folderWritable(const char* path);

    bool fReadable(std::string &f);
    bool fReadable(const char * f);
    void rptprintargs(int argc, char** argv, FILE *rptfile);    
#include "TTimer.hh"
#include <errno.h>
#include <cstring>
#include <sys/stat.h>
    
    
    inline size_t nwrite(int fd, const void *buf, size_t count)
    {
        if(count == 0)
            return count;
        TTimer t(REAL_TIME);
        t.start();
        float mt = 0.0;
        size_t i = 0;
        int errNum;
        errno = 0;
        do 
        {
            ssize_t  wb = write(fd, (char *)buf + i, count - i);
            if(wb >= 0)
                i += wb;
            if(i >= count)
                return i;
            errNum = errno;
            if(!(errNum == ENOSPC || errNum == EINTR || errNum == 0))
                break;
            t.stop();
            if(t.diff() > 30*60)   // 30 minutes of "No Space" errors.
                break;
            if(t.diff() > mt + 30) {
                fprintf(stderr,"File IO Error: %s.\n", strerror(errNum));
                mt = t.diff();
            }
            if(errNum == ENOSPC)
                sleep(30);
        } while(true);
        errno = errNum;
        return i;
    }
    
    inline size_t nfwrite (const void *buf, size_t sz, size_t count, FILE *f)
    {
        if(count == 0)
            return count;
        TTimer t(REAL_TIME);
        t.start();
        float mt = 0.0;
        size_t i = 0;
        int errNum = 0;
        do 
        {
            i += fwrite((char *)buf + sz * i, sz, count - i, f);
            if(i >= count)
                return i;
            errNum = errno;
            if(errNum != ENOSPC && errNum != EINTR)
                break;
            t.stop();
            if(t.diff() > 30*60)   // 30 minutes of "No Space" errors.
                break;
            if(t.diff() > mt + 30) {
                fprintf(stderr,"File IO Error: %s.\n", strerror(errNum));
                mt = t.diff();
            }
            if(errNum == ENOSPC)
                sleep(30);
        } while(true);
        errno = errNum;
        return i;
    }

    bool ispipe(int fno);
        
    bool ispipe(std::string &fname);
#endif
    
