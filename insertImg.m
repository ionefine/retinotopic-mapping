function img = insertImg(img,testImg)
% insertImg(img,testImg)
% 
% Inserts testImg into the center of img.  
% if testImg is larger than img, testImg is cropped and centered.

if size(testImg,1)>size(img,1)
    x0 = ceil((size(testImg,1)-size(img,1))/2)+1;
    testImg = testImg(x0:(x0+size(img,1)-1),:);
end

if size(testImg,2)>size(img,2)
    y0 = ceil((size(testImg,2)-size(img,2))/2)+1;
    testImg = testImg(:,y0:(y0+size(img,2)-1),:);
end


x0 = ceil((size(img,2)-size(testImg,2))/2)+1;
y0 = ceil((size(img,1)-size(testImg,1))/2)+1;
img(y0:(y0+size(testImg,1)-1),x0:(x0+size(testImg,2)-1)) = testImg;