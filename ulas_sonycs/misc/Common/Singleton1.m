classdef Singleton1 < handle
    % Abstract class for the design pattern SINGLETON
    
     methods(Abstract, Static)      
      % If the object doesn't exist or is invalid, create it otherwise return the existing one 
      % in persistent memory.      
      obj = instance();
   end
    
end

