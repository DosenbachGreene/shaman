classdef NullDataProvider < DataProvider
    % A DataProvider that does not provide any data.
    %
    % Useful as a default value wherever the abstract DataProvider class
    % must be instantiated.
    methods
        function data = dataAt(~, ~)
            error("The null data provider cannot provide any data.");
        end
       
        function [data, index] = nextData(~)
            error("The null data provider cannot provide any data.");
        end
        
        function moreData = isMoreData(~)
            moreData = false;
        end
        
        function reset(~)
        end
    end
end