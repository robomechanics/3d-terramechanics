% Copyright 2022 Robomechanics Lab, Carnegie Mellon University
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function colorOutput = cmuColor(printType)
% Returns a 3x1 vector of RBG values for the requested official CMU color 
% given in 'printType'. No argument will return 'red-web'. The currently 
% supported colors areprimary colors 'red-web', 'red-print', 'gray', and 
% 'dark-gray', and secondary colors 'gold', 'teal', 'blue', 'green', and 
% 'dark-green'. 'random' will select one of these colors at random. Nathan
% added 'carnegie-red','black','steel-gray','iron-gray','scots-rose','gold-thread','green-thread',
% 'teal-thread','blue-thread','sky-blue','hall-tan','brick-beige',
% 'hornbostel-teal','palladian-green','weaver-blue','skibo-red'

if nargin == 0
    printType = 'web';
end

colorList = {'red-web', 'red-print', 'gray', 'dark-gray', 'gold',...
    'teal', 'blue', 'green', 'dark-green','carnegie-red','black',...
    'steel-gray','iron-gray','scots-rose','gold-thread','green-thread',...
    'teal-thread','blue-thread','sky-blue','hall-tan','brick-beige',...
    'hornbostel-teal','palladian-green','weaver-blue','skibo-red'};

if strcmp(printType, 'random')
    printType = colorList{randi(length(colorList))};
end

switch printType
    case 'red-web'
        colorOutput = 1/255*[187 0 0];
    case 'red-print'
        colorOutput = 1/255*[176 28 46];
    case 'gray'
        colorOutput = 1/255*[244 244 244];
    case 'dark-gray'
        colorOutput = 1/255*[102 102 102];
    case 'gold'
%         colorOutput = 1/255*[170 102 0];
colorOutput = 1/255*[255 173 0];
    case 'teal'
        colorOutput = 1/255*[0 102 119];
    case 'blue'
        colorOutput = 1/255*[34 68 119];
    case 'green'
        colorOutput = 1/255*[18 220 0];
%         colorOutput = 1/255*[0 136 85];
    case 'dark-green'
        colorOutput = 1/255*[34 68 51];
    %%% Nathan's colors
    % Main colors
    case 'carnegie-red'
        cmyk_val = [0 100 79 20];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[166 25 46];
    case 'black'
        cmyk_val = [0 0 0 100];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[0 0 0];
    case 'steel-gray'
        cmyk_val = [0 0 0 30];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[99 102 106];
    case 'iron-gray'
        cmyk_val = [0 0 0 70];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[187 188 188];
    % Secondary colors
    case 'scots-rose'
        cmyk_val = [0 92 72 0];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[239 51 64];
    case 'gold-thread'
        cmyk_val = [0 32 100 0];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[242 160 0];
    case 'green-thread'
        cmyk_val = [92 2 100 12];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[0 132 61];
    case 'teal-thread'
        cmyk_val = [100 0 40 20];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[0 125 138];
    case 'blue-thread'
        cmyk_val = [100 80 6 32];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[0 45 114];
    case 'sky-blue'
        cmyk_val = [100 10 3 16];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[0 130 186];
    % campus palette
    case 'hall-tan'
        cmyk_val = [15 15 30 15];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[183 176 156];
    case 'brick-beige'
        cmyk_val = [10 11 23 0];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[209 204 189];
    case 'hornbostel-teal'
        cmyk_val = [85 50 58 41];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[13 82 87];
    case 'palladian-green'
        cmyk_val = [60 25 45 0];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[120 159 144];
    case 'weaver-blue'
        cmyk_val = [97 84 44 40];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[0 43 73];
    case 'skibo-red'
        cmyk_val = [15 100 87 35];
%         colorOutput = cmyk(cmyk_val);
        colorOutput = 1/255*[138 43 43];
end

 
        function rgb_val = cmyk(cmyk_val)
            cc = cmyk_val(1)/100; mm = cmyk_val(2)/100; 
            yy = cmyk_val(3)/100; kk = cmyk_val(4)/100; 
            R = (1-cc)*(1-kk);
            G = (1-mm)*(1-kk);
            B = (1-yy)*(1-kk);
            rgb_val = [R,G,B];
        end
end