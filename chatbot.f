      program test
! created by   R. Papesch 
! created on   02 Nov 1998 or earlier
! license      GNU GPLV2.0, please see 
!              https://github.com/papesch/fortran/blob/master/LICENSE
! description  Artificial Intelligence program to beat the Turing test      
      
      character*4 reply
      integer age
      print *,'Hello there.'
      read *,reply
      print *, reply,'to you too.  Would you like to talk?'
      read *,reply
      if(reply.eq.'y'.or.reply.eq.'Y') then
         print *,'Sorry I can''t talk, but I can execute my program.'
      else
         print *,'It''s OK, let''s type messages instead.'
      endif
      write(*,5)
 5    format(' I guess this is a pathetic attempt at artificial'
     +,' intelligence.')
      print *,'So, how are you today?',
      read *,reply
      print *,'I guess it''s easy to feel ',reply,' on a day like this'
      print *,'So, how old are you?'
      read *,age
      print *,age,' is way too old for me. I am only 15 minutes old.'
      print *,'Well, that''s it for now.'
      end
