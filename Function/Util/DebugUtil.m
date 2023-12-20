classdef DebugUtil
    properties(Constant)
        EXP = 4;
        WARNING = 3;
        INFO = 2;
        DEBUG = 1;
    end
    properties
        debugLevel;
        fileID;
        filePath;
    end

    methods
        function obj = DebugUtil(debugLevel, out)
            obj.debugLevel = debugLevel;
            obj.filePath = out;
            if strcmp(out, "terminal")
                obj.fileID = 1;
            else
                obj.fileID = fopen(out, 'w');
                fclose(obj.fileID);
            end
            
        end

        function debug(obj, prefix, message)
            if obj.fileID ~= 1
                obj.fileID = fopen(obj.filePath, 'a');
            end
            if(obj.debugLevel <= obj.DEBUG)
                if isstring(message)
                    fprintf(obj.fileID, prefix + "DEBUG:" + message+ "\n");
                elseif iscell(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message{row, coll}));
                            fprintf(obj.fileID, ", ");
                        end
                        fprintf(obj.fileID, "\n");
                    end
                elseif isnumeric(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message(row, coll)));
                        end
                        fprintf(obj.fileID, "\n");
                    end
                end
            end
            if obj.fileID ~= 1
                fclose(obj.fileID);
            end
        end

        function info(obj, prefix, message)
            if obj.fileID ~= 1
                obj.fileID = fopen(obj.filePath, 'a');
            end
            if(obj.debugLevel <= obj.INFO)
                if isstring(message)
                    fprintf(obj.fileID, prefix + "INFO:" + message+ "\n");
                elseif iscell(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message{row, coll}));
                        end
                        fprintf(obj.fileID, "\n");
                    end
                elseif isnumeric(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message(row, coll)));
                        end
                        fprintf(obj.fileID, "\n");
                    end
                end
            end
            if obj.fileID ~= 1
                fclose(obj.fileID);
            end
        end

        function warning(obj, prefix, message)
            if obj.fileID ~= 1
                obj.fileID = fopen(obj.filePath, 'a');
            end
            if(obj.debugLevel <= obj.WARNING)
                if isstring(message)
                    fprintf(obj.fileID, prefix + "WARNING:" + message+ "\n");
                elseif iscell(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message{row, coll}));
                        end
                        fprintf(obj.fileID, "\n");
                    end
                elseif isnumeric(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message(row, coll)));
                        end
                        fprintf(obj.fileID, "\n");
                    end
                end
            end
            if obj.fileID ~= 1
                fclose(obj.fileID);
            end
        end

        function exp(obj, prefix, message)
            if obj.fileID ~= 1
                obj.fileID = fopen(obj.filePath, 'a');
            end
            if(obj.debugLevel <= obj.EXP)
                if isstring(message)
                    fprintf(obj.fileID, prefix + "EXP:" + message+ "\n");
                elseif iscell(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message{row, coll}));
                        end
                        fprintf(obj.fileID, "\n");
                    end
                elseif isnumeric(message)
                    for row = 1:size(message, 1)
                        fprintf(obj.fileID, prefix);
                        for coll = 1:size(message, 2)
                            fprintf(obj.fileID, "%s\t", string(message(row, coll)));
                        end
                        fprintf(obj.fileID, "\n");
                    end
                end
            end
            if obj.fileID ~= 1
                fclose(obj.fileID);
            end
        end
    end
end