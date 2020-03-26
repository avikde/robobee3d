using Images, FileIO

function darkGhost(img)
	imgs = cat(img...; dims=3)
	imgmin = minimum(imgs, dims=3)
	save("test.jpg", imgmin[:,:,1])
end

function darkGhostDir(dirname)
	img = [load(string(dirname,i,".jpg")) for i=1:3]
	imgs = cat(img...; dims=3)
	imgmin = minimum(imgs, dims=3)
	save("test.jpg", imgmin[:,:,1])
end

# darkGhost([
# 	load("../../../Desktop/mod1 4b1/200v 155hz 1.jpg"),
# 	load("../../../Desktop/mod1 4b1/200v 155hz 2.jpg"),
# 	load("../../../Desktop/mod1 4b1/200v 155hz 3.jpg")
# ])
# darkGhostDir("../../../Desktop/hb1 a1/")
# darkGhostDir("../../../Desktop/hb1 4b1/")
# darkGhostDir("../../../Desktop/mod1a1  a/")
# darkGhostDir("../../../Desktop/bborig/")
darkGhostDir("../../../Desktop/bb4l/")
