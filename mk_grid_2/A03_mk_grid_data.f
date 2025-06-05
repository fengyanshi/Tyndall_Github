         parameter(m=257,n=171)
         real x(m,n),y(m,n),cori(m,n),dep(m,n),dep_swan(m,n)

          m1=m-1
          n1=n-3
         
          open(1,file='xx_curv.txt')
           do j=1,n1
            read(1,*)(x(i,j),i=1,m1)
           enddo
          close(1)

          open(1,file='yy_curv.txt')
           do j=1,n1
            read(1,*)(y(i,j),i=1,m1)
           enddo
          close(1)

          open(1,file='dep_curv_correct.txt')
           do j=1,n1
            read(1,*)(dep(i,j),i=1,m1)
           enddo
          close(1)

         depmax=100.
         do j=1,n1
           do i=1,m1
             dep(i,j)=-dep(i,j)
             if(dep(i,j).gt.100)dep(i,j)=100.0
           enddo
         enddo


         depmin=0.1
         do j=1,n1
           do i=1,m1
             dep_swan(i,j)=dep(i,j)
             cori(i,j)=30.1
           enddo
         enddo

! for swan
       do j=1,n1
        do i=1,m1
          if(dep_swan(i,j).lt.0.1) dep_swan(i,j)=0.1
        enddo
        enddo

        open(1,file='x_circ.txt')
         do j=1,n1
           write(1,100)(x(i,j),i=1,m1)
         enddo
        close(1)

        open(1,file='y_circ.txt')
         do j=1,n1
           write(1,100)(y(i,j),i=1,m1)
         enddo
        close(1)

        open(1,file='cori.txt')
         do j=1,n1
           write(1,100)(cori(i,j),i=1,m1)
         enddo
        close(1)

        open(1,file='dep_circ.txt')
         do j=1,n1
           write(1,100)(dep(i,j),i=1,m1)
         enddo
        close(1)

! for swan
         open(1,file='xxyy_swan.txt')
         do j=1,n1
         write(1,100)(x(i,j),i=1,m1)
         enddo
         do j=1,n1
         write(1,100)(y(i,j),i=1,m1)
         enddo
         close(1)

         open(1,file='dep_swan.txt')
         do j=1,n1
         write(1,100)(dep_swan(i,j),i=1,m1)
         enddo
         close(1)

100     format(3000f12.3)

        end
