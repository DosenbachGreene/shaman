classdef NullDataProvider < DataProvider
    % A DataProvider that does not provide any data.
    %
    % Useful as a default value for properties of the abstract DataProvider
    % class that must be instantiated.
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