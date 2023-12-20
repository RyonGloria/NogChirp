classdef PyramidDecoder
    properties
        loraSet;
        writeDir;
        CollisionFilePath;
        matlabFlagFilePath;
        grloraFlagfilePath;
        binFIlePath;
        binRecord;
    end

    methods
        function obj = PyramidDecoder(loraSet)
            obj.loraSet = loraSet;
            obj.writeDir = "E:\Pyramid_samples\conflict\";
            obj.CollisionFilePath = obj.writeDir + "collision.sigmf-data";
            obj.matlabFlagFilePath = obj.writeDir + "matlabFlag.txt";
            obj.grloraFlagfilePath = obj.writeDir + "grloraFlag.txt";
            obj.binFIlePath = obj.writeDir + "bin.txt";
        end

        function obj = decode(obj, signals)
            obj.binRecord = cell(1, 0);
            % 删除目录下所有文件
            obj.clearDir(obj.writeDir);
            % 写入读取区
            obj.writeSignalToFile(signals, obj.CollisionFilePath);
            % 写入标志位
            writematrix([], obj.matlabFlagFilePath);
            % 等待gr-lora处理完成
            % fprintf("wait gr-lora...\n");
            while ~exist(obj.grloraFlagfilePath, 'file')
                pause(1);
            end
            % fprintf("gr-lora Done!\n");
            % 清除gr-lora标志文件
            delete(obj.grloraFlagfilePath);
            try
                obj = obj.readBin();
            catch
                obj.binRecord = cell(1, 0);
            end
        end

        function obj = clearDir(obj, dirPath)
            % 删除文件夹下所有内容
            if exist(dirPath, 'dir')
                % 获取文件夹中的所有文件
                files = dir(fullfile(dirPath, '*'));
                % 删除文件夹中的所有文件
                for i = 1:length(files)
                    if ~files(i).isdir
                        filePath = fullfile(dirPath, files(i).name);
                        delete(filePath);
                    end
                end
            else
                disp(['文件夹 ' dirPath ' 不存在。']);
            end
        end

        function obj = writeSignalToFile(obj, signals, fileName)
            signalLength = length(signals)*2;
            signalProcessed = zeros(1, signalLength);
            signalReal = real(signals);   signalImag = imag(signals);
            signalProcessed(1:2:signalLength-1) = signalReal;
            signalProcessed(2:2:signalLength) = signalImag;

            fid=fopen(fileName, 'wb');
            fwrite(fid, signalProcessed, 'float32');
            fclose(fid);
        end

        function obj = readBin(obj)
            % 打开文件
            fileID = fopen(obj.binFIlePath, 'r');
            % 使用textscan读取每行内容为字符串
            data_cell_array = textscan(fileID, '%s', 'Delimiter', '\n');
            % 关闭文件
            fclose(fileID);
            % 获取每行的字符串
            data_strings = data_cell_array{1};
            % 初始化一个元胞数组来保存数据
            obj.binRecord = cell(numel(data_strings), 1);
            % 将每行的字符串分割并放入元胞数组中
            for i = 1:numel(data_strings)
                % 使用strsplit将逗号分隔的字符串拆分为单元格数组
                line_data = strsplit(data_strings{i}, ',');
            
                % 将空字符串替换为-1，并将该行数据存储在元胞数组的相应位置，并转换为double数组
                line_data(cellfun(@isempty, line_data)) = {'-1'};
                obj.binRecord{i} = cellfun(@str2double, line_data) + 1;
            end
        end
    end
end