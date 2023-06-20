function data1 =     insertstimframeback(data1, data2,stim_frames,insertionmethod)
   if strcmp(insertionmethod, 'meanofprepost')
       for i = 1:length(stim_frames)
           if stim_frames(i) > size(data1,2)
              data1 = cat(2, data1, data2(:,end));
           elseif stim_frames(i) == 1
              data1 = cat(2, data1(:,1), data1);
           else
              data1 = cat(2, data1(:,1:stim_frames(i)-1), mean(data1(:,stim_frames(i)-1:stim_frames(i)),2),  data1(:,stim_frames(i):end));
           end
       end
   else
       
        for i = 1:length(stim_frames)
           if stim_frames(i) > size(data1,2)
              data1 = cat(2, data1, data2(:,stim_frames(i)));
           elseif stim_frames(i) == 1
              data1 = cat(2, data2(:,stim_frames(i)), data1);

           else
              data1 = cat(2, data1(:,1:stim_frames(i)-1), data2(:,stim_frames(i)),  data1(:,stim_frames(i):end));
           end
        end
   end
end

