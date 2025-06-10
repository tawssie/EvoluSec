function JArray=CharVector2JavaArray(CharVector)

JArray=javaArray('java.lang.Character',numel(CharVector));

for i=1:numel(CharVector)
    JArray(i)=java.lang.Character(CharVector(i));
end

